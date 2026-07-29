#!/usr/bin/env python3
"""Publish MF radar quick-look plots and host telemetry to juha.no."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import tempfile
import time


PROCESS_PATTERNS = (
    ("ringbuffer", "Digital RF ring buffer", "drf ringbuffer"),
    ("recorder_10", "USRP recorder · 192.168.10.2", "thor.py -m 192.168.10.2"),
    ("recorder_20", "USRP recorder · 192.168.20.2", "thor.py -m 192.168.20.2"),
    ("processor", "Cross-correlation processor", "python calc_xc.py"),
)
CHANNELS = ("ch1", "ch2", "ch3", "ch4")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="/data1/mfraw")
    parser.add_argument("--metadata-dir", default="/data2/metadata/xc")
    parser.add_argument("--plot-dir", default="/data2/plots/monitor")
    parser.add_argument("--remote", default="j@juha.no:/var/www/html/mf/")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def iso_utc(timestamp: float) -> str:
    return (
        dt.datetime.fromtimestamp(timestamp, tz=dt.timezone.utc)
        .isoformat()
        .replace("+00:00", "Z")
    )


def read_meminfo() -> dict:
    values = {}
    with Path("/proc/meminfo").open(encoding="utf-8") as handle:
        for line in handle:
            key, raw = line.split(":", 1)
            values[key] = int(raw.strip().split()[0]) * 1024
    total = values["MemTotal"]
    available = values["MemAvailable"]
    return {
        "total_bytes": total,
        "available_bytes": available,
        "used_percent": 100.0 * (total - available) / total,
    }


def process_rows() -> list[dict]:
    processes = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            command = (entry / "cmdline").read_bytes().replace(b"\0", b" ").decode().strip()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if command:
            processes.append((int(entry.name), command))

    result = []
    for name, label, pattern in PROCESS_PATTERNS:
        matches = [(pid, command) for pid, command in processes if pattern in command]
        result.append(
            {
                "name": name,
                "label": label,
                "alive": bool(matches),
                "pid": matches[0][0] if matches else None,
            }
        )
    timer_active = subprocess.run(
        ["systemctl", "is-active", "--quiet", "mf-radar-monitor.timer"],
        check=False,
    ).returncode == 0
    result.append(
        {
            "name": "monitor_timer",
            "label": "Web monitor refresh timer",
            "alive": timer_active,
            "pid": None,
        }
    )
    wind_timer_active = subprocess.run(
        ["systemctl", "is-active", "--quiet", "mf-radar-winds.timer"],
        check=False,
    ).returncode == 0
    result.append(
        {
            "name": "wind_timer",
            "label": "Automatic wind processor timer",
            "alive": wind_timer_active,
            "pid": None,
        }
    )
    return result


def disk_rows() -> list[dict]:
    roles = {
        "/": "Operating system",
        "/data1": "Raw-voltage ring buffer",
        "/data2": "Analyzed wind products",
        "/data3": "Available data disk",
        "/data4": "Available data disk",
    }
    rows = []
    for mount, role in roles.items():
        usage = shutil.disk_usage(mount)
        rows.append(
            {
                "mount": mount,
                "role": role,
                "total_bytes": usage.total,
                "free_bytes": usage.free,
                "used_percent": 100.0 * usage.used / usage.total,
            }
        )
    return rows


def newest_file(directory: Path, pattern: str) -> Path | None:
    candidates = []
    for path in directory.glob(pattern):
        try:
            candidates.append((path.stat().st_mtime, path))
        except FileNotFoundError:
            # Digital RF files can be atomically renamed while telemetry is read.
            continue
    return max(candidates)[1] if candidates else None


def acquisition_status(raw_dir: Path) -> tuple[dict, float, int]:
    channels = {}
    ages = []
    now = time.time()
    for channel in CHANNELS:
        channel_dir = raw_dir / channel
        date_dirs = sorted(path for path in channel_dir.glob("20*") if path.is_dir())
        latest = newest_file(date_dirs[-1], "rf@*.h5") if date_dirs else None
        age = max(0.0, now - latest.stat().st_mtime) if latest else 1e12
        channels[channel] = {"acquisition_age_seconds": age}
        ages.append(age)
    return channels, max(ages), sum(age < 30 for age in ages)


def processing_lag(metadata_dir: Path) -> float:
    date_dirs = sorted(path for path in metadata_dir.glob("20*") if path.is_dir())
    latest = newest_file(date_dirs[-1], "xc@*.h5") if date_dirs else None
    if latest is None:
        return 1e12
    match = re.search(r"xc@(\d+)", latest.name)
    return max(0.0, time.time() - int(match.group(1))) if match else 1e12


def read_rti_status(plot_dir: Path) -> dict | None:
    path = plot_dir / "rti_status.json"
    if not path.exists():
        return None
    value = json.loads(path.read_text(encoding="utf-8"))
    window_end = dt.datetime.fromisoformat(value["window_end"].replace("Z", "+00:00"))
    value["age_seconds"] = max(0.0, time.time() - window_end.timestamp())
    return value


def build_status(args: argparse.Namespace) -> dict:
    now = time.time()
    channels, acquisition_age, channels_ok = acquisition_status(Path(args.raw_dir))
    processes = process_rows()
    lag = processing_lag(Path(args.metadata_dir))
    rti = read_rti_status(Path(args.plot_dir))
    memory = read_meminfo()
    disks = disk_rows()

    if not all(process["alive"] for process in processes) or channels_ok < 4:
        overall = "error"
    elif lag > 300 or rti is None or rti["age_seconds"] > 600 or memory["used_percent"] > 90:
        overall = "warning"
    else:
        overall = "ok"

    try:
        uptime_seconds = float(Path("/proc/uptime").read_text().split()[0])
    except (OSError, ValueError):
        uptime_seconds = 0.0

    return {
        "schema_version": 1,
        "generated_at": iso_utc(now),
        "host": socket.gethostname(),
        "overall_status": overall,
        "uptime_seconds": uptime_seconds,
        "acquisition_age_seconds": acquisition_age,
        "acquisition_channels_ok": channels_ok,
        "processing_lag_seconds": lag,
        "rti": rti,
        "memory": memory,
        "disks": disks,
        "processes": processes,
        "channels": channels,
    }


def main() -> None:
    args = parse_args()
    plot_dir = Path(args.plot_dir)
    required = (
        plot_dir / "latest_rti_48h_full.png",
        plot_dir / "latest_rti_48h_mesosphere.png",
    )
    optional = (
        plot_dir / "latest_selected_wind_pixels.png",
        plot_dir / "latest_altitude_cuts_30m.png",
        plot_dir / "latest_snr_doppler_30m_0_200.png",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Monitor plots are missing: {', '.join(missing)}")

    with tempfile.TemporaryDirectory(prefix="mf-radar-monitor-") as temporary:
        staging = Path(temporary)
        for source in required:
            shutil.copy2(source, staging / source.name)
        for source in optional:
            if source.exists():
                shutil.copy2(source, staging / source.name)
        (staging / "status.json").write_text(
            json.dumps(build_status(args), indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        command = [
            "rsync",
            "-az",
            "--delay-updates",
            "--chmod=D755,F644",
            f"{staging}/",
            args.remote,
        ]
        if args.dry_run:
            command.insert(1, "--dry-run")
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()

const REFRESH_MS = 30_000;
const plots = [
  {
    id: "rti_full",
    title: "Ionospheric propagation",
    label: "Full range · 0–1500 km",
    file: "latest_rti_48h_full.png",
  },
  {
    id: "rti_mesosphere",
    title: "Partial reflections & meteor echoes",
    label: "Mesosphere · 0–200 km",
    file: "latest_rti_48h_mesosphere.png",
  },
  {
    id: "rti_selected",
    title: "Pixels accepted for wind processing",
    label: "Automatic selection · 50–200 km",
    file: "latest_selected_wind_pixels.png",
  },
];

const dialog = document.querySelector("#plot-dialog");
let activePlotIndex = 0;
let latestStatus = null;

function formatDuration(seconds) {
  if (!Number.isFinite(seconds)) return "—";
  if (seconds < 60) return `${Math.max(0, Math.round(seconds))} s`;
  if (seconds < 3600) return `${Math.round(seconds / 60)} min`;
  if (seconds < 86400) return `${(seconds / 3600).toFixed(seconds < 10_800 ? 1 : 0)} h`;
  return `${(seconds / 86400).toFixed(1)} d`;
}

function formatBytes(bytes) {
  if (!Number.isFinite(bytes)) return "—";
  const units = ["B", "KiB", "MiB", "GiB", "TiB"];
  let value = bytes;
  let unit = 0;
  while (value >= 1024 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(unit > 2 ? 1 : 0)} ${units[unit]}`;
}

function formatUtc(iso) {
  const date = new Date(iso);
  if (Number.isNaN(date.valueOf())) return "Unknown";
  return new Intl.DateTimeFormat("en-GB", {
    dateStyle: "medium",
    timeStyle: "medium",
    timeZone: "UTC",
  }).format(date) + " UTC";
}

function stateLabel(state) {
  return {
    ok: "Operational",
    warning: "Degraded",
    error: "Offline",
    unknown: "Connecting",
  }[state] || "Unknown";
}

function updateClock() {
  const text = new Intl.DateTimeFormat("en-GB", {
    hour: "2-digit",
    minute: "2-digit",
    second: "2-digit",
    hour12: false,
    timeZone: "UTC",
  }).format(new Date());
  document.querySelector("#utc-clock").textContent = `${text} UTC`;
}

function renderProcesses(processes = []) {
  const target = document.querySelector("#process-list");
  target.replaceChildren();
  processes.forEach((process) => {
    const row = document.createElement("div");
    row.className = "process-row";
    row.innerHTML = `
      <span class="process-indicator ${process.alive ? "alive" : ""}" aria-hidden="true"></span>
      <span class="process-name">${process.label}</span>
      <span class="process-meta">${process.alive ? (process.pid ? `PID ${process.pid}` : "active") : "not running"}</span>
    `;
    target.append(row);
  });
}

function renderStorage(disks = []) {
  const target = document.querySelector("#storage-list");
  target.replaceChildren();
  disks.forEach((disk) => {
    const usedPercent = Math.round(disk.used_percent);
    const row = document.createElement("div");
    row.className = "storage-row";
    row.innerHTML = `
      <span class="storage-name">${disk.role || disk.mount} <small>${disk.mount}</small></span>
      <span class="storage-meta">${usedPercent}% · ${formatBytes(disk.free_bytes)} free</span>
      <div class="storage-bar" aria-label="${disk.mount} ${usedPercent}% used">
        <span style="width: ${Math.min(100, usedPercent)}%"></span>
      </div>
    `;
    target.append(row);
  });
}

function renderPlot(plot, status, cacheKey) {
  const rti = status.rti;
  if (!rti) return;

  const age = rti.age_seconds;
  const freshness = document.querySelector(`#${plot.id}-freshness`);
  freshness.textContent = age < 120 ? "Current" : `${formatDuration(age)} old`;
  freshness.classList.toggle("warning", age > 300);
  document.querySelector(`#${plot.id}-time`).textContent =
    `Window ends · ${formatUtc(rti.window_end)}`;

  const image = document.querySelector(`#${plot.id}-image`);
  image.src = `${plot.file}?v=${cacheKey}`;
}

function render(status) {
  latestStatus = status;
  const state = status.overall_status || "unknown";
  const pill = document.querySelector("#overall-pill");
  pill.dataset.state = state;
  document.querySelector("#overall-label").textContent = stateLabel(state);
  document.querySelector("#last-update").textContent =
    `Telemetry ${formatDuration(status.telemetry_age_seconds)} ago`;

  const acquisitionAge = status.acquisition_age_seconds;
  document.querySelector("#acquisition-value").textContent =
    acquisitionAge < 30 ? "Streaming" : `${formatDuration(acquisitionAge)} old`;
  document.querySelector("#acquisition-detail").textContent =
    `${status.acquisition_channels_ok}/4 channels current`;

  const lag = status.processing_lag_seconds;
  document.querySelector("#lag-value").textContent = formatDuration(lag);
  document.querySelector("#lag-detail").textContent =
    lag > 300 ? "Processor is behind ring-buffer head" : "Processing near realtime";

  const memory = status.memory;
  document.querySelector("#memory-value").textContent =
    `${Math.round(memory.used_percent)}% used`;
  document.querySelector("#memory-detail").textContent =
    `${formatBytes(memory.available_bytes)} available`;

  const data1 = status.disks.find((disk) => disk.mount === "/data1");
  if (data1) {
    document.querySelector("#disk-value").textContent =
      `${Math.round(data1.used_percent)}% used`;
    document.querySelector("#disk-detail").textContent =
      `${formatBytes(data1.free_bytes)} available`;
  }

  document.querySelector("#host-value").textContent = status.host;
  document.querySelector("#uptime-value").textContent = formatDuration(status.uptime_seconds);

  const cacheKey = Date.parse(status.generated_at) || Date.now();
  plots.forEach((plot) => renderPlot(plot, status, cacheKey));
  renderProcesses(status.processes);
  renderStorage(status.disks);

  if (dialog.open) updateDialog(activePlotIndex);
}

async function refresh() {
  try {
    const response = await fetch(`status.json?v=${Date.now()}`, { cache: "no-store" });
    if (!response.ok) throw new Error(`Status request failed: ${response.status}`);
    const status = await response.json();
    status.telemetry_age_seconds = Math.max(
      0,
      (Date.now() - Date.parse(status.generated_at)) / 1000,
    );
    render(status);
  } catch (error) {
    const pill = document.querySelector("#overall-pill");
    pill.dataset.state = "error";
    document.querySelector("#overall-label").textContent = "No telemetry";
    document.querySelector("#last-update").textContent = "Status feed unavailable";
    console.error(error);
  }
}

function updateDialog(index) {
  activePlotIndex = (index + plots.length) % plots.length;
  const plot = plots[activePlotIndex];
  const sourceImage = document.querySelector(`#${plot.id}-image`);
  document.querySelector("#dialog-channel").textContent = plot.label;
  document.querySelector("#dialog-title").textContent = plot.title;
  document.querySelector("#dialog-image").src = sourceImage.src;
  const dataTime = latestStatus?.rti?.window_end;
  document.querySelector("#dialog-time").textContent =
    dataTime ? formatUtc(dataTime) : "Awaiting plot metadata";
}

function openDialog(index) {
  updateDialog(index);
  dialog.showModal();
}

document.querySelectorAll("[data-plot-index]").forEach((button) => {
  button.addEventListener("click", () => openDialog(Number(button.dataset.plotIndex)));
});

document.querySelector("#dialog-close").addEventListener("click", () => dialog.close());
document.querySelector("#dialog-prev").addEventListener("click", () => updateDialog(activePlotIndex - 1));
document.querySelector("#dialog-next").addEventListener("click", () => updateDialog(activePlotIndex + 1));

dialog.addEventListener("click", (event) => {
  if (event.target === dialog) dialog.close();
});

document.addEventListener("keydown", (event) => {
  if (event.key === "Escape" && dialog.open) {
    dialog.close();
    return;
  }
  if (event.key === "ArrowLeft") {
    if (dialog.open) updateDialog(activePlotIndex - 1);
    else document.querySelectorAll("[data-plot-index]")[0].focus();
  }
  if (event.key === "ArrowRight") {
    if (dialog.open) updateDialog(activePlotIndex + 1);
    else document.querySelectorAll("[data-plot-index]")[1].focus();
  }
});

updateClock();
setInterval(updateClock, 1000);
refresh();
setInterval(refresh, REFRESH_MS);

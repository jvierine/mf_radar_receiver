import numpy as n
import datetime
import time
from datetime import datetime, timedelta, timezone

xc_dir="/data2/metadata/xc"
raw_dir="/data1/mfraw"
dc_offset=-0.25-0.25j

def unix2date(x):
    # ATC bug-fix for python3
    return datetime.datetime.utcfromtimestamp(float(x))

def unix2dirname(x):
    return (unix2date(x).strftime('%Y-%m-%d'))

def fir_lowpass_hann(fc=20e3, fs=1000000, num_taps=100):
    # Normalized cutoff (0 to 1)
    norm_cutoff = fc / (fs / 2)

    # Ideal sinc filter
    x = n.arange(num_taps) - (num_taps - 1) / 2
    h = n.sinc(norm_cutoff * x)

    # Hann window
    w = n.hanning(num_taps)
    h = h * w

    # Normalize
    h = h / n.sum(h)
    return h
def unix2datestr(unix_seconds):
    date_string = datetime.datetime.utcfromtimestamp(unix_seconds).isoformat() + "Z"
    return(date_string)


def get_prev_day_bounds():
    # Get current time in UTC
    now = datetime.now(timezone.utc)

    # Previous day
    yesterday = now.date() - timedelta(days=1)
    
    # Start of previous day (00:00:00 UTC)
    start = datetime(yesterday.year, yesterday.month, yesterday.day, tzinfo=timezone.utc)
    
    # End of previous day (23:59:59 UTC)
    end = start + timedelta(days=1) - timedelta(seconds=1)
    
    start_unix = int(start.timestamp())
    end_unix = int(end.timestamp())
    return(start_unix,end_unix)


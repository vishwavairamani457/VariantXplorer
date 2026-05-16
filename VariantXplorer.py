import os
import sys
import time
import threading
import subprocess
import requests
import webview

def start_streamlit():
    """Start Streamlit using the bundled Python executable"""
    subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "frontend.py",
         "--server.headless=true",
         "--server.port=8501"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

def wait_for_streamlit():
    """Keep checking until Streamlit is ready to serve"""
    print("Waiting for Streamlit to start...")
    while True:
        try:
            response = requests.get("http://localhost:8501")
            if response.status_code == 200:
                print("Streamlit is ready!")
                break
        except Exception:
            pass
        time.sleep(0.5)  # Retry every 0.5 seconds

if __name__ == '__main__':
    # Step 1: Start Streamlit in background thread
    streamlit_thread = threading.Thread(target=start_streamlit, daemon=True)
    streamlit_thread.start()

    # Step 2: Wait until Streamlit server is fully ready
    wait_for_streamlit()

    # Step 3: Open the app in a native window
    window = webview.create_window(
        title="VariantXplorer",
        url="http://localhost:8501",
        width=1280,
        height=800,
        resizable=True,
        min_size=(800, 600)
    )

    webview.start()
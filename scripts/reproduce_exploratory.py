#!/usr/bin/env python3
"""Run the independent exploratory RSV and climate analysis."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from extensions.humidity.exploratory.runner import main

if __name__ == "__main__":
    main()

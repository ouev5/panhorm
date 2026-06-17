#!/bin/bash
cd /www/wwwroot/multi-omics
source venv/bin/activate
export PYTHONPATH=/www/wwwroot/multi-omics/src:$PYTHONPATH
python -m src.webui.app

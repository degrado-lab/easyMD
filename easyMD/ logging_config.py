# my_md_package/logging_config.py

import logging
import sys

logger = logging.getLogger("easyMD")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[%(levelname)s] %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
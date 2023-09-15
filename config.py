import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_FOLDER = os.path.join(BASE_DIR, "input")
#INPUT_FOLDER = os.environ.get("INPUT_FOLDER", INPUT_FOLDER)

OUTPUT_FOLDER = os.path.join(BASE_DIR, "output")
#OUTPUT_FOLDER = os.environ.get("OUTPUT_FOLDER", OUTPUT_FOLDER)
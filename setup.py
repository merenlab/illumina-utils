import glob
from setuptools import setup


def get_scripts():
    # Preserve the previous behavior: all scripts in scripts/, except those ending with '-OBSOLETE'
    return [
        script for script in glob.glob("scripts/*")
        if not script.endswith("-OBSOLETE")
    ]


if __name__ == "__main__":
    setup(scripts=get_scripts())

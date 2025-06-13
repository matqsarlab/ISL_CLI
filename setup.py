from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="ISL",
    version="0.2.0",
    py_modules=["isl"],
    packages=find_packages(include=["lib", "alvaDescCLIWrapper"]),
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "isl=isl:main",
        ],
    },
    python_requires=">=3.12",
)

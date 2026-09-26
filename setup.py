from setuptools import find_packages, setup


setup(
    name="step-upgrade",
    version="0.3.0",
    description="Scalable storm identification and lineage tracking for precipitation grids.",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=["numpy>=1.18,<3.0", "scipy>=1.4,<2.0"],
    extras_require={"envelope": ["scikit-image>=0.19,<1.0"]},
)

from setuptools import find_packages, setup


setup(
    name="step-upgrade",
    version="0.1.0",
    description="Scalable storm identification and tracking for precipitation grids.",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=["numpy>=1.18", "scipy>=1.4", "scikit-image>=0.17"],
)

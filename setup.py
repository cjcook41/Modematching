import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

setuptools.setup(
	name = "modematch",
	version = "0.0.7",
	author="Cameron Cook",
	author_email="cam.cook32@gmail.com",
	description="Crystal free energy prediction suite",
	long_discription=long_description,
	long_description_content_type="text/markdown",
	packages=setuptools.find_packages(),
	install_requires=[
		'numpy',
		'sympy',
		'pandas',
		'matching',
		],
	url="https://github.com/cjcook41/Modematching.git",
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
	],
	python_requires='>=3.7',
)

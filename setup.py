from setuptools import setup

def readme():
	with open('README.md', 'r') as f:
		return f.read()

setup(name='extremevariantfilter',
      version='0.0a2',
      description='A set of tools to aid in the identification of false positive variants in Variant Call Files.',
      long_description=readme(),
      long_description_content_type="text/markdown",
      classifiers=[
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Intended Audience :: Science/Research',
         'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      url='https://github.com/Ellis-Anderson/extremevariantfilter',
      author='Complete Genomics',
      author_email='eanderson@genomics.cn',
      packages=['extremevariantfilter'],
      install_requires=[
         'scikit-learn',
         'numpy',
         'pandas',
         'xgboost',
         'docopt'
      ],
      scripts=[
         'bin/train_model',
         'bin/apply_filter'],
      include_package_data=True)

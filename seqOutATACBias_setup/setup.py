import setuptools

long_description = """
seqOutATACBias first generates a union bedGraph of all included
seqOutBias masks, using the "masks" command, then applies a rule
ensemble model to correct Tn5 bias using the "scale" command.

Please consult the included README.md file for installation and
dependency information.
"""

setuptools.setup(
    name='seqOutATACBias',
    version='1.0',
    scripts=['./scripts/seqOutATACBias',
            './scripts/maskgeneration',
            './scripts/RE_implementation.R'],
    author='Jacob Wolpe',
    description='seqOutATACBias: Rule Ensemble Modeling of seqOutBias Scaling',
    packages=['seqOutATACBias'],
    install_requires=[
        'setuptools'
    ],
    python_requires='>=3.9.12'
)

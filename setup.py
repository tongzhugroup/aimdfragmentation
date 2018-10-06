from setuptools import setup
setup(name='AIMDFragmentation',
      version='1.0.6',
      description='AIMD Fragmentation Calculation',
      keywords="AIMD Fragmentation",
      url='https://github.com/njzjz/AIMDFragmentation',
      author='Jinzhe Zeng',
      author_email='jzzeng@stu.ecnu.edu.cn',
      packages=['AIMDFragmentation','AIMDBlock'],
      install_requires=['numpy','GaussianRunner'])

from setuptools import setup

with open("README.rst", "r") as fh: 
    description = fh.read() 
    print(description)
setup(
    name='KonfersiPy',
    version='0.1.0',    
    description='Konfersi Indonesia Python Library',
    url='https://github.com/Konfersi-Indonesia/KonfersiPy.git',
    author='Konfersi Indonesia',
    author_email='konfersi.indonesia@gmail.com',
    license='MIT License',
    packages=['KonfersiPy'],
    install_requires=['numpy','matplotlib','xarray','basemap','pandas','rasterio','netcdf4','fiona','shapely','geopandas']
)

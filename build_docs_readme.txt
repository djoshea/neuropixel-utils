conda create --name mkdocs pip
conda activate mkdocs
pip install mkdocs mkdocs-material

mkdocs serve # test local
mkdocs gh-deploy # deploy to master


find -name '*.py' -exec black -l 80 '{}' ';'
find -name '*.py' -exec flake8 --ignore=E501,E241,W503 '{}' ';'


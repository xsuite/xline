find -name '*.py' -exec black '{}' ';'
find -name '*.py' -exec flake8 --ignore=E501,E241,W503 '{}' ';'


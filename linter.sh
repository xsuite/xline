find -name '*.py' -exec black -l 79 '{}' ';'
flake8 --ignore=E501,E241,W503
pytest --cov=pysixtrack --cov-report html

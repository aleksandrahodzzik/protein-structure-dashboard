IMAGE ?= protein-structure-dashboard
PORT ?= 8501

.PHONY: run test lint format requirements docker-build docker-run

run:  ## start the app locally on http://localhost:8501
	uv run streamlit run streamlit_app.py

test:  ## unit tests with coverage (no network needed)
	uv run pytest --cov --cov-report=term-missing

lint:  ## ruff lint + format check + mypy, exactly what CI runs
	uv run ruff check .
	uv run ruff format --check .
	uv run mypy

format:  ## auto-fix lint issues and reformat
	uv run ruff check --fix .
	uv run ruff format .

requirements:  ## regenerate requirements.txt for Streamlit Community Cloud
	uv export --frozen --no-dev --no-hashes --no-emit-project -o requirements.txt
	echo "-e ." >> requirements.txt

docker-build:
	docker build -t $(IMAGE) .

docker-run:  ## run the container, then open http://localhost:$(PORT)
	docker run --rm -p $(PORT):8501 $(IMAGE)

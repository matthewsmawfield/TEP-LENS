.PHONY: install pipeline test manuscript pdf clean help

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install: ## Install Python dependencies
	pip install -r requirements.txt

pipeline: ## Run the full analysis pipeline (steps 00–35)
	python scripts/steps/run_all_steps.py

test: ## Run regression tests
	python -m pytest tests/ -v

manuscript: ## Build manuscript from HTML components
	cd site && npm ci && npm run build

pdf: ## Generate PDF (requires playwright + chromium)
	python scripts/generate_site_pdf.py --quality high --wait-time 5

clean: ## Remove generated outputs and logs
	rm -rf logs/*.log
	@echo "Cleaned logs. Results preserved in results/outputs/."

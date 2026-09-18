# Multi-Crop Contextual Bayesian Optimization Benchmark

This repository tests Contextual Bayesian Optimization (CBO) strategies across 5 different crops (Lettuce, Strawberry, Spinach, Asparagus, Cabbage). It directly compares Ivan's Heteroskedastic CBO, Standard CBO, Standard BO, and Random Search.

### Benchmark Setup
- **Strict Budget:** All models are capped at exactly 20 physical trials to test their extreme data efficiency.
- **Cost-Aware:** Standard and Contextual BO utilize an Expected Loss risk filter dynamically tied to the remaining budget.
- **Spatial Diversity:** Uses a "Search-then-Relax" EIC batching mechanism.

### Running via GitHub Actions
This repository is configured to run fully in the cloud via GitHub Actions to save local CPU cycles:
1. Push to your repository.
2. Go to the **Actions** tab.
3. Click on **Multi-Crop BO Benchmark** and hit **Run workflow**.
4. Download the `multi-crop-results` artifact zip when finished to view the CSVs for all 5 crops.

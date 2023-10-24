
# Order of Execution (Files Only)

1. `ProteoCorrAnalyser.py`
2. `BayesianWeightOptimization.py`
3. `ProximityBasedCysFeSelector.py`
4. `ANMModeAveragesCompiler.py`
5. `island_dynamic_programming_V9.py`
6. `TileBasedSimilarityVisualization.R`
7. `ProteinTileSimViz.R`
8. `ProteinChainOptimalPathAnalyzer.R`

## Example execution:

1\. **`ProteoCorrAnalyser.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProteoCorrAnalyser.py 5t5i 2fdn All G A;
   ```

2\. **`BayesianWeightOptimization.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/BayesianWeightOptimization.py All 5t5i.chainG.nowa
   ```

3\. **`ProximityBasedCysFeSelector.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProximityBasedCysFeSelector.py All 5t5i G;
   ```

4\. **`ANMModeAveragesCompiler.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ANMModeAveragesCompiler.py All
   ```

5\. **`island_dynamic_programming_V9.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/island_dynamic_programming_V9.py
   ```

6\. **`TileBasedSimilarityVisualization.R`**

   ```bash
   Rscript --vanilla /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/TileBasedSimilarityVisualization.R 5t5i "Chain G" All
   ```

7\. **`ProteinTileSimViz.R`**

   ```bash
   Rscript --vanilla /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProteinTileSimViz.R 5t5i "Chain G" All
   ```

8\. **`ProteinChainOptimalPathAnalyzer.R`**

   ```bash
   Rscript --vanilla /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProteinChainOptimalPathAnalyzer.R --pdb-id=5t5i --chain-id=G
   ```


Explanation:
- The `#` symbol is used for headings.
- The `1.`, `2.`, `3.`, etc. are used for ordered lists.
- The backticks (`` ` ``) are used for inline code formatting.
- Triple backticks (```bash ... ```) are used for multi-line code blocks.
- Bold text is created using double asterisks (`**`).
# MDS-FSM: Coverage-Based Frequent Subgraph Mining in Single Graphs

A high-performance C++20 implementation of the MDS-FSM algorithm for mining top-k frequent subgraph patterns using the Minimum Density Support (MDS) measure.

## Requirements

- CMake 3.20 or later
- A C++20 compiler
- AVX2 and POPCNT support are recommended for best performance

## Usage

```bash
# Build first
cd MDS-FSM
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# Run the binary
./build/mds-fsm -d <dataset> -k <k> [-o <output_dir>]
```

## Parameters

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `-d <dataset>` | Path to the dataset file (.gfu format), or dataset name in `dataset/` folder | Yes | - |
| `-k <k>` | Number of top-k patterns to mine | Yes | - |
| `-o <output_dir>` | Output directory for results | No | `output/` |

## Examples

```bash
# Mine top-10 patterns from yeast dataset (using dataset name)
./build/mds-fsm -d yeast -k 10

# Mine top-10 patterns from yeast dataset (using full path)
./build/mds-fsm -d dataset/yeast.gfu -k 10

# Mine top-20 patterns with custom output directory
./build/mds-fsm -d dataset/yeast.gfu -k 20 -o results/
```

## Input Format

The input file uses the GFU (Graph File Unified) format:

```text
t # <graph_id>
v <vertex_id> <vertex_label>
...
e <source_id> <target_id> <edge_label>
...
```

## Output Format

The program outputs:

1. Dataset statistics (vertices, edges, labels)
2. Total runtime
3. Top-k patterns with their MDS values

Each pattern is represented as a DFS code sequence:

```text
( from, to, from_label, edge_label, to_label )
```

Results are also saved to a text file in the output directory.

## Datasets

The `dataset/` directory contains example dataset:

- `yeast.gfu` - Yeast protein interaction network

## Example

```bash
$ ./build/mds-fsm -d yeast -k 10

Dataset: dataset/yeast.gfu
Vertices: 3112, Edges: 12519
Unique vertex labels: 71, Unique edge labels: 1

Total runtime: 0.032s
==================
      Answer
==================
Size of answer set: 10
Minimum MDS: 155
--------
Pattern 1: MDS = 155
( 0, 1, 2, 1, 2 )
( 1, 2, 2, 1, 36 )
--------
Pattern 2: MDS = 158
( 0, 1, 2, 1, 36 )
--------
Pattern 3: MDS = 162
( 0, 1, 0, 1, 2 )
( 1, 2, 2, 1, 2 )
( 2, 3, 2, 1, 2 )
--------
Pattern 4: MDS = 164
( 0, 1, 2, 1, 2 )
( 1, 2, 2, 1, 2 )
( 2, 3, 2, 1, 3 )
--------
Pattern 5: MDS = 171
( 0, 1, 2, 1, 2 )
( 1, 2, 2, 1, 3 )
...
```

## Citation

If you use this code in your research, please cite:

```bibtex
@inproceedings{mdsfsm2026,
  title={MDS-FSM: Coverage-Based Frequent Subgraph Mining in Single Graphs},
  author={Guo, Xiaozhen and Liu, Xueli and Dong, Bowen and Wan, Li and Ge, Jiake and Ma, Shuai},
  booktitle={Proceedings of the VLDB Endowment},
  year={2026}
}
```

## License

This project is released for academic research purposes.

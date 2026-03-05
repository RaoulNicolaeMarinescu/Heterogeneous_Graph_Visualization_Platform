# Heterogeneous Graph Visualization Platform

This project was developed as my Bachelor's thesis in Digital Communication at the University of Milan.

The goal of the project is to design and implement an interactive platform for the visualization and exploration of heterogeneous biological graphs representing relationships between genes, Gene Ontology (GO) terms, and Human Phenotype Ontology (HPO) terms.

The system focuses on usability and human–computer interaction principles to make complex biological data easier to explore and interpret.

---

# Overview

Modern biological datasets often involve complex relationships between genes, biological functions, and phenotypes. These relationships are typically represented through ontologies such as:

- Gene Ontology (GO)
- Human Phenotype Ontology (HPO)

However, these datasets are difficult to explore using traditional tabular representations.

This project proposes a visual exploration system that allows users to interactively navigate these relationships through a graph-based interface.

The platform allows users to:

- explore relationships between genes and biological functions
- navigate ontology hierarchies
- visualize connections between genes, GO terms, and HPO terms
- progressively expand the graph to reduce cognitive overload

The project focuses on **data visualization, usability, and interactive exploration of complex datasets** rather than biological analysis itself.

The system was developed at **AnacletoLab – Department of Computer Science, University of Milan**. :contentReference[oaicite:0]{index=0}

---

# Problem

Biological datasets describing gene functions and phenotypes can involve thousands of nodes and hierarchical relationships.

For example:

- genes are associated with GO terms describing biological functions
- genes are associated with HPO terms describing phenotypic effects
- GO and HPO terms themselves form complex hierarchical ontologies

These structures are difficult to interpret because:

- ontologies contain tens of thousands of nodes
- hierarchical relationships span many levels
- traditional representations (matrices, CSV files) are not visually intuitive

The goal of the project is therefore **not to discover new biological insights**, but to create a tool that helps users explore these relationships interactively and intuitively. :contentReference[oaicite:1]{index=1}

---

# System Architecture

The platform follows a clear separation between two main components:

### Backend (Python)

Responsible for:

- preprocessing biological datasets
- parsing ontology files
- constructing heterogeneous graphs
- exporting graph structures as JSON files

The backend processes multiple data sources including:

- Gene–GO matrices
- Gene–HPO matrices
- Gene Ontology (.obo)
- Human Phenotype Ontology (.obo)

The output is a set of JSON files describing nodes and edges of the graph.

### Frontend (Web Application)

The frontend is an interactive web application built using:

- HTML
- CSS
- JavaScript
- D3.js

It is responsible for:

- rendering the graph visualization
- handling user interaction
- dynamically expanding the graph based on user input
- managing different exploration modes

The interface supports zoom, pan, node selection, and progressive expansion of relationships. :contentReference[oaicite:2]{index=2}

---

# Data Pipeline

The backend performs several preprocessing steps:

### 1. Parsing ontology files

Ontology files in `.obo` format are parsed to extract:

- node identifiers
- labels
- hierarchical relationships (`is_a`, `part_of`)

This information is used to reconstruct the ontology hierarchy.

### 2. Building the heterogeneous graph

The system generates nodes representing:

- Genes
- GO terms
- HPO terms

Edges represent different types of relationships:

- GO → GO (hierarchical relationships)
- HPO → HPO (hierarchical relationships)
- Gene → GO
- Gene → HPO

### 3. Graph export

The final graph is exported as JSON:

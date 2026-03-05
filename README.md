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

nodes: { id, label, type, depth }
edges: { source, target, relation }

These JSON files are then loaded by the frontend visualization.

---

# Visualization Approach

Displaying the entire graph simultaneously would produce a visually overwhelming interface.

To address this problem, the system adopts **progressive disclosure**, meaning that the graph expands only when users interact with specific nodes.

Three main exploration modes are implemented:

### Gene-centric mode

Users start from a gene node and explore:

- associated GO terms
- associated HPO terms
- ontology hierarchies

### GO-centric mode

Users start from a GO term and explore:

- hierarchical relationships within the ontology
- associated genes

### HPO-centric mode

Users start from a phenotype term and explore:

- phenotype hierarchy
- related genes
- associated biological functions

This interaction model allows users to gradually explore the graph without visual overload. :contentReference[oaicite:3]{index=3}

---

# Interface Design

The interface was designed following usability principles from human–computer interaction research.

Key design choices include:

- progressive disclosure of information
- hierarchical node layout
- color-coded node types
- zoom and pan navigation
- contextual information panels
- light/dark theme support

Nodes are organized spatially according to their type and hierarchical depth to maintain visual clarity. :contentReference[oaicite:4]{index=4}

---

# Technologies Used

Backend

- Python
- Pandas
- CSV processing
- JSON generation

Frontend

- HTML
- CSS
- JavaScript
- D3.js

Data Sources

- Gene Ontology
- Human Phenotype Ontology
- gene–GO association matrices
- gene–HPO association matrices

---

# Features

- Interactive graph visualization
- Multiple exploration modes
- Ontology hierarchy reconstruction
- Progressive graph expansion
- Dynamic filtering of nodes
- Node information tooltips
- Light and dark interface modes

---

# Limitations

The project focuses on data exploration rather than biological analysis.

Limitations include:

- no statistical or predictive analysis
- scalability limitations for very large datasets
- no integration with external biological databases

These aspects represent potential directions for future improvements.

---

# Future Work

Possible future developments include:

- user testing and usability evaluation
- integration with external biological data sources
- support for larger datasets
- personalized exploration views
- saving and sharing exploration sessions

---

# Author

Raoul Nicolae Marinescu  
Bachelor’s Degree in Digital Communication  
University of Milan

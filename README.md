# micron-analysis-challange

Methodology
1. Convert Demand to Tool Requirements

For each quarter:

Take wafer loading for Node1, Node2, and Node3.

Multiply loading by required process time per step (RPT).

Sum total processing minutes required per workstation type.

Convert required minutes into number of tools using:

Available minutes per week

Tool utilization rates

Round up to ensure sufficient capacity.

This guarantees that demand can be processed within the given utilization limits.

2. Search for Feasible Production Splits

The solver explores different production splits:

Node1 split between Fab1 and Fab2

Node2 split between Fab1 and Fab2

Node3:

C-steps can be assigned to any fab

Non-C steps are split between Fab1 and Fab2

For each candidate allocation:

Compute tools required per fab

Calculate space used per fab

Compare against fab space limits:

Fab1: 1500 m²

Fab2: 1300 m²

Fab3: 700 m²

3. Select the Best Valid Plan

Among all feasible allocations:

Choose the configuration with the lowest total space usage

If no solution satisfies space constraints, return a fallback configuration

This ensures manufacturability while prioritizing space efficiency.

4. Outputs

For each quarter, the model:

Prints production allocation per fab

Reports tool counts per workstation per fab

Displays space usage versus capacity

Exports:

A detailed flow plan (per node, per step, per fab)

Exact tool counts per workstation per fab

Additional Validation

The script also:

Computes total required space if production were consolidated

Checks if total required space fits within the combined 3-fab capacity (3500 m²)

Evaluates high-growth future quarters for feasibility

Summary

This is a capacity-driven allocation model that:

Demand → Processing Minutes → Tool Count → Space Usage → Feasibility Check

It systematically searches for a feasible and space-efficient production plan under strict physical constraints.

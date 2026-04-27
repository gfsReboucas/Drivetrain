# Drivetrain Modeling and Scaling Framework

This project implements a physics-based framework for modeling and scaling multi-stage geared drivetrains while preserving their structural and dynamic characteristics (e.g., frequency distribution).

## Why this matters
Engineering models are often either:
- Too detailed (computationally expensive), or
- Too simplified (losing physical accuracy)

This project addresses that gap by enabling reduced-order drivetrain representations that remain physically meaningful.

## What it does
- Models spur and helical gear systems (parallel and planetary)
- Preserves dynamic behavior during drivetrain scaling
- Enables integration with industry tools:
  - [Simpack](https://www.3ds.com/products-services/simulia/products/simpack/) (multibody dynamics)
  - [KISSsoft](https://www.kisssoft.com/en) (gear design and safety factors)

## Engineering value
- Supports faster system-level simulations
- Enables early-stage design exploration
- Bridges detailed gear design and system-level dynamics

## Tech
MATLAB, object-oriented design, Jupyter notebooks, COM interface integration

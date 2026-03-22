# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository processes historical EOSS (Edge of Space Sciences) high-altitude balloon flight data. It queries APRS telemetry packets from a PostgreSQL database, computes flight physics (velocity, acceleration, air density, Reynolds number transitions), and produces multiple output formats.

## Running

### Rust (primary)

```bash
cargo build --release
cargo run --release                                    # all flights, all outputs
cargo run --release -- --flight EOSS-391               # single flight
cargo run --release -- --flight EOSS-391 --output-type parquet  # single flight, parquet only
cargo run --release -- --dbname eosstracker            # use production DB
```

Requires a local PostgreSQL database with the `packets` table (schema in `legacy-database.sql`). The database must have PostGIS extensions (`geometry` types). Database name defaults to `legacy` (configurable via `--dbname` or `EOSS_DBNAME` env var). Uses `rayon` for parallel flight processing.

### Python (legacy, kept for reference)

```bash
python3 process-data.py
```

### Python Dependencies (legacy only)

pandas, numpy, scipy, matplotlib, psycopg2, pytz, simplekml, xlsxwriter

## Architecture

### Rust Implementation

Multi-module design in `src/`:
- **main.rs** — CLI entry (clap), rayon orchestration, output dispatch
- **config.rs** — CLI args: `--dbname`, `--flight`, `--output-type`, `--flightlist`, `--output-dir`
- **models.rs** — Serde structs for FlightInput (from JSON) and FlightMetadata (enriched output)
- **db.rs** — PostgreSQL query via `postgres` crate, DST detection via `chrono-tz`, DataFrame construction
- **processing.rs** — Core pipeline: per-beacon processing, consolidation, trimming, ascent/descent split, Reynolds detection, polynomial fitting
- **physics.rs** — Haversine distance, air density, polynomial fitting (nalgebra SVD), VMR degree selection
- **output/** — CSV (Polars CsvWriter), JSON (serde_json), Parquet (Polars ParquetWriter, replaces pickle), XLSX (rust_xlsxwriter), KML (quick-xml), metadata consolidation

### Python Implementation (legacy)

**Single-script design** — all logic lives in `process-data.py`.

### Data Flow

1. `main()` reads `flightlist.json` (flight metadata: beacons, dates, weights, parachute info)
2. For each flight, spawns a `multiprocessing.Process` running `processThread()`
3. `processThread()` calls:
   - `queryDatabase()` — SQL query against the `legacy` PostgreSQL database, extracts APRS packets with parsed telemetry (temperature, pressure from raw packet strings)
   - `process_df()` — core analysis: trims pre-launch/post-landing data using moving averages, splits into ascent/descent, computes velocity/acceleration/air density, detects Reynolds number transitions via distance-to-line analysis, fits polynomial curves (degree auto-selected via VMR)
   - `createPlot()` — 4-panel matplotlib figure (altitude vs time, velocity vs altitude, acceleration, ACF)
   - `createKML()` — Google Earth KML with 3D flight paths, waypoints at 10k ft intervals, and Reynolds transition markers
4. After all processes complete, consolidates per-flight JSON into `output/json/flights_metadata.json` and summary CSV/XLSX

### Output Directory Structure

All outputs go to `output/` subdirectories: `csv/`, `json/`, `parquet/` (Rust) or `pkl/` (Python), `png/` (Python only), `xlsx/`, `kml/`

### Key Domain Concepts

- **Reynolds transitions**: detected by finding points furthest from the average ascent line (distance_to_line), indicating shifts between laminar (low Re) and turbulent (high Re) airflow regimes
- **Flight trimming**: pre-launch packets removed using short (3-pt) and long (10-pt) forward moving averages of vertical velocity; post-landing stragglers removed by time gap detection
- **Outlier ejection**: ascent rates constrained to (-5, 50) ft/s, descent rates to < 5 ft/s
- **Polynomial curve fitting**: degree selected adaptively based on Variance-to-Mean Ratio (VMR) of vertical velocity, max degree 13
- **Timezone handling**: `is_dst()` determines MST vs MDT for SQL query time windows (Mountain time zone)

### Database

PostgreSQL with PostGIS. The `packets` table stores raw APRS packets with geometry columns for 2D/3D locations. Temperature and pressure are parsed from raw packet strings using SQL regex patterns in the query.

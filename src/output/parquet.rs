use anyhow::Result;
use polars::prelude::*;
use std::fs::File;

/// Write DataFrame to Parquet format (replaces pickle)
pub fn write_parquet(df: &DataFrame, path: &str) -> Result<()> {
    let mut file = File::create(path)?;
    ParquetWriter::new(&mut file).finish(&mut df.clone())?;
    Ok(())
}

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_permanova::{Config, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-permanova", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Square distance-matrix TSV (lsmat layout: blank top-left, labeled rows/cols).
    input: PathBuf,

    /// Grouping table: `id<tab>group`, one sample per line (header optional).
    #[arg(short = 'g', long)]
    grouping: PathBuf,

    /// Permutations for the p-value; 0 skips it (p-value = NA).
    #[arg(short = 'p', long, default_value_t = 999)]
    permutations: usize,

    /// Treat inputs as comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Decimal places in the output.
    #[arg(long, default_value_t = 6)]
    precision: usize,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let delim = if self.csv { ',' } else { '\t' };
        let cfg = Config {
            permutations: self.permutations,
            seed: self.common.seed.unwrap_or(0),
            threads: self.common.thread_count(),
            delim,
            precision: self.precision,
        };

        let dm_reader =
            BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?);
        let grouping_reader = BufReader::new(File::open(&self.grouping).map_err(|e| {
            RsomicsError::InvalidInput(format!("{}: {e}", self.grouping.display()))
        })?);
        let mut out: Box<dyn Write> = if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };
        run(dm_reader, grouping_reader, &mut out, &cfg)?;
        out.flush().map_err(RsomicsError::Io)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "PERMANOVA pseudo-F test for group differences from a distance matrix.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.stats.distance.permanova",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1111/j.1442-9993.2001.01070.pp.x"),
    }),
    usage_lines: &["dm.tsv -g groups.tsv [-p 999] [--seed S] [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('g'),
                long: "grouping",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Grouping table: id<tab>group per line (header optional).",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "permutations",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("999"),
                description: "Permutations for the p-value (0 skips it).",
                why_default: Some("scikit-bio default; precision 1/(1+999)=0.001"),
            },
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse inputs as comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "precision",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("6"),
                description: "Decimal places in the output.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "PERMANOVA with 999 permutations and a fixed seed",
            command: "rsomics-permanova dm.tsv -g groups.tsv --seed 42",
        },
        Example {
            description: "Statistic only (no permutation p-value)",
            command: "rsomics-permanova dm.tsv -g groups.tsv -p 0",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}

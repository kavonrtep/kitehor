use thiserror::Error;

#[derive(Debug, Error)]
pub enum HordetectError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTA parse error: {0}")]
    Fasta(String),

    #[error("invalid parameter: {0}")]
    InvalidParam(String),

    #[error("array too short: {id} has {length} bp, need at least {min} bp")]
    ArrayTooShort {
        id: String,
        length: usize,
        min: usize,
    },

    #[error("array {id} is {n_fraction:.1}% N (limit: {limit:.1}%)")]
    TooManyNs {
        id: String,
        n_fraction: f64,
        limit: f64,
    },

    #[error("{0}")]
    Other(String),
}

pub type Result<T> = std::result::Result<T, HordetectError>;

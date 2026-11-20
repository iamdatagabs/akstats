"""Command-line interface."""
import argparse
import logging
from akstats.pipeline import run_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the AKStats exploratory pipeline.")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )
    run_pipeline()


if __name__ == "__main__":
    main()

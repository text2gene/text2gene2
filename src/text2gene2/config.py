from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="ignore")

    # NCBI
    ncbi_api_key: str = ""

    # Google CSE
    google_api_key: str = ""
    google_cse_whitelist_id: str = "009155410218757639293:nyh3tdzvfhc"
    google_cse_schema_id: str = "003914143621252222636:-mop04_esug"

    # Redis
    redis_url: str = "redis://localhost:6379/0"

    # Cache TTLs (seconds)
    cache_ttl_lvg: int = 604800
    cache_ttl_litvar2: int = 86400
    cache_ttl_clinvar: int = 86400
    cache_ttl_google: int = 604800
    cache_ttl_europepmc: int = 86400

    # VariantValidator
    variant_validator_url: str = "https://rest.variantvalidator.org"
    variant_validator_genome_build: str = "GRCh38"
    variant_validator_rate_limit: float = 1.5  # req/sec

    # App
    app_host: str = "0.0.0.0"
    app_port: int = 8000
    app_debug: bool = False


settings = Settings()

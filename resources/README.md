# Resources

Raw input data for the pipeline — everything the workflow reads *before* it
computes anything. Two subfolders by provenance:

| Subfolder | Contents | See |
|-----------|----------|-----|
| [`public/`](public/) | Third-party datasets from open sources (UN Comtrade trade records, World Bank unit-value indices). Pre-seeded so the whole pipeline runs without any API key. | [`public/README.md`](public/README.md) |
| [`inhouse/`](inhouse/) | Data curated by the authors for this project — the FAO ↔ HS product correspondence table. | [`inhouse/README.md`](inhouse/README.md) |

> The pre-seeded files in `public/` let every downstream rule run out of the box.
> Only the `get_uncomtrade` rule needs the `comtrade_apikey` environment variable,
> and only if you want to re-download the raw trade data.

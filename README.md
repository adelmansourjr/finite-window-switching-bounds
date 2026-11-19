# Finite Window Switching Bounds

Finite-window joint-occupancy bounds for a simulated two-switch system, using known mode fields and an envelope for the switching implementation. The archived evaluation covered 1,032 trajectories. The report compares interval width with product and maximum-entropy forecasts and gives the limitations.

[Report](SWITCHING_REPORT.md) · [Protocol](SWITCHING_PROTOCOL.md) · [Archived results](switching_results.json)

Run with the Python 3.10+ standard library.

```sh
python3 switching.py --verify switching_results.json
```

This public release extracts the original experiment from the Bristol research workspace. It contains generated trajectories and aggregate outcomes. [Related telemetry study](https://github.com/adelmansourjr/switching-telemetry).

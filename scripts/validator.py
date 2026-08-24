#!/usr/bin/env python3
"""
scripts/validator.py

Validates Vina docking energy predictions against biophysical reference constants.
Checks if predicted binding energy shifts (ΔΔG = ΔG_T790M - ΔG_WT) correctly separate
1st, 2nd, and 3rd generation EGFR inhibitors based on resistance profiles.
"""

# THIS IS MEANT TO FAIL, SINCE CHANGE IN CHANGE IN BINDING ENERGY ISNT THE BEST WAY TO SHOW
#, PLIP DOES A BETTER JOB, INLCUDE THIS IN THE PAPER

import json
import sys
from pathlib import Path


def load_json(file_path: Path) -> dict:
    if not file_path.exists():
        print(f"[-] Error: Required file not found at {file_path}")
        sys.exit(1)
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_affinity(drug_data: dict, target_key: str) -> float | None:
    """Extracts affinity values across common key naming variations."""
    possible_keys = [
        f"mean_delta_g_{target_key}",
        f"{target_key}_delta_g",
        f"{target_key}_affinity",
        f"delta_g_{target_key}",
        f"mean_{target_key}",
        target_key
    ]
    for k in possible_keys:
        if k in drug_data and isinstance(drug_data[k], (int, float)):
            return float(drug_data[k])

    if target_key in drug_data and isinstance(drug_data[target_key], dict):
        sub_dict = drug_data[target_key]
        for sub_k in ["mean_affinity", "affinity", "delta_g", "mean_delta_g", "mean"]:
            if sub_k in sub_dict and isinstance(sub_dict[sub_k], (int, float)):
                return float(sub_dict[sub_k])

    return None


def validate_docking_results():
    eval_file = Path("data/processed/docking_results/docking_evaluation.json")
    if not eval_file.exists():
        eval_file = Path("results/docking_evaluation.json")
    if not eval_file.exists():
        eval_file = Path("docking_evaluation.json")

    benchmark_file = Path("results/biophysical_reference_benchmark.json")

    eval_data = load_json(eval_file)
    benchmark_data = load_json(benchmark_file)

    THRESHOLDS = {
        "1st_gen_failure_min": 0.3,
        "2nd_gen_retention_max": 0.3,
        "2nd_gen_retention_min": 0.1,
        "3rd_gen_selectivity_max": 0.1,
    }

    print("\n" + "=" * 65)
    print("      EGFR INHIBITOR DOCKING VALIDATION REPORT")
    print("=" * 65 + "\n")

    summary_results = []

    # Handle both top-level drug dictionary and nested "drugs" / "results" key
    drugs = eval_data.get("drugs", eval_data.get("results", eval_data))
    ref_drugs = benchmark_data.get("drugs", benchmark_data)

    parsed_drug_count = 0

    for drug_name, metrics in drugs.items():
        if not isinstance(metrics, dict):
            continue

        drug_key = drug_name.lower()
        if drug_key not in ref_drugs:
            print(f"[!] Warning: {drug_name} not found in reference benchmark. Skipping.")
            continue

        ref_info = ref_drugs[drug_key]
        expected_gen = ref_info.get("generation", "Unknown")

        # Flexible extraction for WT and T790M values
        delta_g_wt = extract_affinity(metrics, "wt")
        if delta_g_wt is None:
            delta_g_wt = extract_affinity(metrics, "wildtype")

        delta_g_mut = extract_affinity(metrics, "t790m")
        if delta_g_mut is None:
            delta_g_mut = extract_affinity(metrics, "mutant")

        ddg_explicit = metrics.get("ddg", metrics.get("delta_delta_g", metrics.get("mutation_shift")))

        if delta_g_wt is not None and delta_g_mut is not None:
            ddg = delta_g_mut - delta_g_wt
        elif ddg_explicit is not None:
            ddg = float(ddg_explicit)
        else:
            print(f"[-] Error: Could not extract ΔG or ΔΔG values for {drug_name}.")
            continue

        parsed_drug_count += 1

        # Classify drug behavior based on predicted ΔΔG
        if ddg >= THRESHOLDS["1st_gen_failure_min"]:
            predicted_gen = "1st Gen"
            status_note = "Severe affinity loss against T790M (Resistance)"
        elif THRESHOLDS["2nd_gen_retention_min"] <= ddg < THRESHOLDS["2nd_gen_retention_max"]:
            predicted_gen = "2nd Gen"
            status_note = "Moderate affinity shift"
        else:  # ddg < 0.1 kcal/mol
            predicted_gen = "3rd Gen"
            status_note = "Sustained/Enhanced binding against T790M"

        is_correct = (predicted_gen == expected_gen)
        status_flag = "PASS" if is_correct else "FAIL"

        ref_ki = ref_info.get("ki_nm", {})

        summary_results.append({
            "drug": drug_key,
            "expected_gen": expected_gen,
            "predicted_gen": predicted_gen,
            "ddg": round(ddg, 2),
            "wt_delta_g": round(delta_g_wt, 2) if delta_g_wt is not None else None,
            "mut_delta_g": round(delta_g_mut, 2) if delta_g_mut is not None else None,
            "status": status_flag,
            "note": status_note,
            "ref_ki": ref_ki
        })

        print(f"Drug: {drug_key.upper()} ({expected_gen})")
        if delta_g_wt is not None and delta_g_mut is not None:
            print(f"  * WT ΔG: {delta_g_wt:.2f} kcal/mol | T790M ΔG: {delta_g_mut:.2f} kcal/mol")
        print(f"  * Calculated ΔΔG: {ddg:+.2f} kcal/mol")
        print(f"  * Behavioral Classification: {predicted_gen} ({status_note})")
        print(f"  * Reference Ki Parameters: {ref_ki}")
        print(f"  * Validation Result: [{status_flag}]\n" + "-" * 65)

    # Determine validation status based on parsed drug count and pass flags
    all_passed = all(item["status"] == "PASS" for item in summary_results) if summary_results else False
    validation_passed = (parsed_drug_count > 0) and all_passed

    output_data = {
        "validation_passed": validation_passed,
        "thresholds_kcal_mol": THRESHOLDS,
        "results": summary_results
    }

    output_dir = Path("results")
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "validation_report.json"

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=4)

    print(f"\n[+] Validation report saved to -> {out_path.resolve()}")

    if validation_passed:
        print(f"[=== SUCCESS ===] Validated {parsed_drug_count} drugs successfully!\n")
        sys.exit(0)
    else:
        print(
            f"[=== FAILURE ===] Processed {parsed_drug_count} drugs. One or more predictions failed classification.\n")
        sys.exit(1)


if __name__ == "__main__":
    validate_docking_results()
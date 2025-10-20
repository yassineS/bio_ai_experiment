#!/usr/bin/env python3
"""
Simple script to generate text-based visualizations of the analysis results.
This creates ASCII charts for the top packages analysis.

Usage:
    python3 visualize_results.py [csv_file]
    
If no csv_file is provided, defaults to 'all_200_packages_ranked.csv'
"""

import csv
import json
import sys
from collections import Counter

# Default data file
DEFAULT_CSV_FILE = 'all_200_packages_ranked.csv'


def load_csv_data(filename):
    """Load package data from CSV file."""
    packages = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            packages.append(row)
    return packages


def create_horizontal_bar(value, max_value, width=40):
    """Create a simple horizontal bar chart."""
    filled = int((float(value) / max_value) * width)
    bar = '█' * filled + '░' * (width - filled)
    return bar


def print_top_packages_chart(packages, top_n=10):
    """Print a chart of top N packages by improvement score."""
    print("\n" + "="*80)
    print(f"TOP {top_n} PACKAGES BY IMPROVEMENT POTENTIAL")
    print("="*80 + "\n")
    
    max_score = float(packages[0]['composite_improvement_score'])
    
    for i, pkg in enumerate(packages[:top_n], 1):
        score = float(pkg['composite_improvement_score'])
        name = pkg['name']
        lang = pkg['language']
        
        bar = create_horizontal_bar(score, max_score, width=50)
        
        print(f"{i:2d}. {name:20s} [{lang:15s}] {bar} {score:.2f}")
    
    print()


def print_language_distribution(packages, top_n=50):
    """Print language distribution for top N packages."""
    print("\n" + "="*80)
    print(f"LANGUAGE DISTRIBUTION (TOP {top_n} PACKAGES)")
    print("="*80 + "\n")
    
    # Extract primary language (before /)
    languages = [pkg['language'].split('/')[0] for pkg in packages[:top_n]]
    lang_counts = Counter(languages)
    
    max_count = max(lang_counts.values())
    
    for lang, count in sorted(lang_counts.items(), key=lambda x: x[1], reverse=True):
        bar = create_horizontal_bar(count, max_count, width=40)
        percentage = (count / top_n) * 100
        print(f"{lang:15s} {bar} {count:2d} packages ({percentage:5.1f}%)")
    
    print()


def print_category_distribution(packages, top_n=50):
    """Print category distribution for top N packages."""
    print("\n" + "="*80)
    print(f"CATEGORY DISTRIBUTION (TOP {top_n} PACKAGES)")
    print("="*80 + "\n")
    
    categories = [pkg['category'] for pkg in packages[:top_n]]
    cat_counts = Counter(categories)
    
    max_count = max(cat_counts.values())
    
    # Sort by count
    for cat, count in sorted(cat_counts.items(), key=lambda x: x[1], reverse=True):
        bar = create_horizontal_bar(count, max_count, width=40)
        percentage = (count / top_n) * 100
        cat_name = cat.replace('_', ' ').title()
        print(f"{cat_name:25s} {bar} {count:2d} ({percentage:5.1f}%)")
    
    print()


def print_quality_metrics(packages, top_n=50):
    """Print average quality metrics for top N packages."""
    print("\n" + "="*80)
    print(f"AVERAGE QUALITY METRICS (TOP {top_n} PACKAGES)")
    print("="*80 + "\n")
    
    top_pkgs = packages[:top_n]
    
    avg_code = sum(float(p['code_quality_score']) for p in top_pkgs) / len(top_pkgs)
    avg_doc = sum(float(p['doc_quality_score']) for p in top_pkgs) / len(top_pkgs)
    avg_test = sum(float(p['test_coverage_estimate']) for p in top_pkgs) / len(top_pkgs)
    avg_pop = sum(float(p['popularity_estimate']) for p in top_pkgs) / len(top_pkgs)
    
    print(f"Code Quality:       {avg_code:5.2f}/10 {create_horizontal_bar(avg_code, 10, 30)} Target: 8.0")
    print(f"Documentation:      {avg_doc:5.2f}/10 {create_horizontal_bar(avg_doc, 10, 30)} Target: 8.0")
    print(f"Test Coverage:      {avg_test:5.2f}/10 {create_horizontal_bar(avg_test, 10, 30)} Target: 8.0")
    print(f"Popularity:         {avg_pop:5.2f}/10 {create_horizontal_bar(avg_pop, 10, 30)}")
    
    print("\nImprovement Gaps:")
    print(f"  Code Quality:     {8.0 - avg_code:+.2f} points needed")
    print(f"  Documentation:    {8.0 - avg_doc:+.2f} points needed")
    print(f"  Test Coverage:    {8.0 - avg_test:+.2f} points needed")
    
    print()


def print_score_distribution(packages):
    """Print distribution of improvement scores."""
    print("\n" + "="*80)
    print("IMPROVEMENT SCORE DISTRIBUTION (ALL PACKAGES)")
    print("="*80 + "\n")
    
    scores = [float(p['composite_improvement_score']) for p in packages]
    
    # Create histogram bins
    bins = [
        (60, 70, "60-70 (Critical)"),
        (55, 60, "55-60 (High)"),
        (50, 55, "50-55 (Medium-High)"),
        (45, 50, "45-50 (Medium)"),
        (40, 45, "40-45 (Low-Medium)"),
        (0, 40, "0-40 (Low)"),
    ]
    
    max_count = 0
    bin_counts = []
    
    for low, high, label in bins:
        count = sum(1 for s in scores if low <= s < high)
        bin_counts.append((label, count))
        max_count = max(max_count, count)
    
    for label, count in bin_counts:
        bar = create_horizontal_bar(count, max_count, width=40)
        percentage = (count / len(packages)) * 100
        print(f"{label:25s} {bar} {count:3d} packages ({percentage:5.1f}%)")
    
    print()


def print_comparison_table(packages):
    """Print comparison of metrics across all packages."""
    print("\n" + "="*80)
    print("QUALITY METRICS COMPARISON")
    print("="*80 + "\n")
    
    scores = {
        'code': [float(p['code_quality_score']) for p in packages],
        'doc': [float(p['doc_quality_score']) for p in packages],
        'test': [float(p['test_coverage_estimate']) for p in packages],
    }
    
    print(f"{'Metric':<20} {'Min':>8} {'Avg':>8} {'Max':>8} {'Median':>8}")
    print("-" * 80)
    
    for name, label in [('code', 'Code Quality'), ('doc', 'Documentation'), ('test', 'Test Coverage')]:
        values = sorted(scores[name])
        min_val = min(values)
        max_val = max(values)
        avg_val = sum(values) / len(values)
        median_val = values[len(values)//2]
        
        print(f"{label:<20} {min_val:8.2f} {avg_val:8.2f} {max_val:8.2f} {median_val:8.2f}")
    
    print()


def main():
    """Main execution function."""
    # Get filename from command line or use default
    csv_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CSV_FILE
    
    # Load data
    print(f"\nLoading package data from {csv_file}...")
    try:
        packages = load_csv_data(csv_file)
        print(f"Loaded {len(packages)} packages\n")
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
        print(f"Usage: python3 {sys.argv[0]} [csv_file]")
        sys.exit(1)
    
    # Generate visualizations
    print_top_packages_chart(packages, top_n=15)
    print_quality_metrics(packages, top_n=50)
    print_language_distribution(packages, top_n=50)
    print_category_distribution(packages, top_n=50)
    print_score_distribution(packages)
    print_comparison_table(packages)
    
    # Summary
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nTotal packages analyzed: {len(packages)}")
    print(f"Top candidates identified: 50")
    print("\nFor detailed information, see:")
    print("  - top_50_packages_for_improvement.md")
    print("  - EXECUTIVE_SUMMARY.md")
    print("  - TOP_10_QUICK_REFERENCE.md")
    print()


if __name__ == "__main__":
    main()

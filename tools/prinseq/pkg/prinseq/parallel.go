package prinseq

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sync"
)

// ParallelResult holds the result of processing a single file
type ParallelResult struct {
	Filename string
	Stats    *Stats
	Error    error
}

// ProcessFilesParallel processes multiple files in parallel
func ProcessFilesParallel(filenames []string, isFastq bool, workers int) ([]ParallelResult, error) {
	if workers <= 0 {
		workers = 4 // Default to 4 workers
	}

	// Create channels
	jobs := make(chan string, len(filenames))
	results := make(chan ParallelResult, len(filenames))

	// Start worker pool
	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for filename := range jobs {
				result := processFile(filename, isFastq)
				results <- result
			}
		}()
	}

	// Send jobs
	for _, filename := range filenames {
		jobs <- filename
	}
	close(jobs)

	// Wait for all workers to finish
	go func() {
		wg.Wait()
		close(results)
	}()

	// Collect results
	allResults := make([]ParallelResult, 0, len(filenames))
	for result := range results {
		allResults = append(allResults, result)
	}

	return allResults, nil
}

// processFile processes a single file
func processFile(filename string, isFastq bool) ParallelResult {
	file, err := os.Open(filename)
	if err != nil {
		return ParallelResult{
			Filename: filename,
			Error:    fmt.Errorf("error opening file: %w", err),
		}
	}
	defer file.Close()

	stats, err := CalculateEnhancedStats(file, isFastq)
	return ParallelResult{
		Filename: filename,
		Stats:    stats,
		Error:    err,
	}
}

// FilterFilesParallel filters multiple files in parallel
func FilterFilesParallel(inputFiles []string, outputDir string, isFastq bool, opts FilterOptions, workers int) error {
	if workers <= 0 {
		workers = 4
	}

	// Create output directory if it doesn't exist
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		return fmt.Errorf("error creating output directory: %w", err)
	}

	// Create channels
	type job struct {
		input  string
		output string
	}
	jobs := make(chan job, len(inputFiles))
	errors := make(chan error, len(inputFiles))

	// Start worker pool
	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				err := filterSingleFile(j.input, j.output, isFastq, opts)
				if err != nil {
					errors <- fmt.Errorf("error processing %s: %w", j.input, err)
				}
			}
		}()
	}

	// Send jobs
	for _, inputFile := range inputFiles {
		outputFile := filepath.Join(outputDir, filepath.Base(inputFile))
		jobs <- job{input: inputFile, output: outputFile}
	}
	close(jobs)

	// Wait for all workers to finish
	wg.Wait()
	close(errors)

	// Check for errors
	var firstError error
	for err := range errors {
		if firstError == nil {
			firstError = err
		}
	}

	return firstError
}

// filterSingleFile filters a single file
func filterSingleFile(inputPath, outputPath string, isFastq bool, opts FilterOptions) error {
	// Open input file
	input, err := os.Open(inputPath)
	if err != nil {
		return fmt.Errorf("error opening input: %w", err)
	}
	defer input.Close()

	// Create output file
	output, err := os.Create(outputPath)
	if err != nil {
		return fmt.Errorf("error creating output: %w", err)
	}
	defer output.Close()

	// Filter
	return Filter(input, output, isFastq, opts)
}

// BatchProcessConfig holds configuration for batch processing
type BatchProcessConfig struct {
	InputFiles  []string
	OutputDir   string
	IsFastq     bool
	Workers     int
	FilterOpts  FilterOptions
	GenerateStats bool
	GenerateReport bool
}

// BatchProcess processes multiple files with various operations
func BatchProcess(config BatchProcessConfig) ([]ParallelResult, error) {
	if config.Workers <= 0 {
		config.Workers = 4
	}

	// Create output directory
	if config.OutputDir != "" {
		if err := os.MkdirAll(config.OutputDir, 0755); err != nil {
			return nil, fmt.Errorf("error creating output directory: %w", err)
		}
	}

	results := make([]ParallelResult, 0, len(config.InputFiles))
	var mu sync.Mutex

	// Create job queue
	type job struct {
		filename string
		index    int
	}
	jobs := make(chan job, len(config.InputFiles))
	
	// Start workers
	var wg sync.WaitGroup
	for i := 0; i < config.Workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				result := processBatchFile(j.filename, config)
				mu.Lock()
				results = append(results, result)
				mu.Unlock()
			}
		}()
	}

	// Queue jobs
	for i, filename := range config.InputFiles {
		jobs <- job{filename: filename, index: i}
	}
	close(jobs)

	// Wait for completion
	wg.Wait()

	return results, nil
}

// processBatchFile processes a single file with all configured operations
func processBatchFile(filename string, config BatchProcessConfig) ParallelResult {
	// Open input file
	input, err := os.Open(filename)
	if err != nil {
		return ParallelResult{
			Filename: filename,
			Error:    fmt.Errorf("error opening file: %w", err),
		}
	}
	defer input.Close()

	// Calculate statistics
	stats, err := CalculateEnhancedStats(input, config.IsFastq)
	if err != nil {
		return ParallelResult{
			Filename: filename,
			Error:    fmt.Errorf("error calculating stats: %w", err),
		}
	}

	result := ParallelResult{
		Filename: filename,
		Stats:    stats,
	}

	// Generate report if requested
	if config.GenerateReport && config.OutputDir != "" {
		reportPath := filepath.Join(config.OutputDir, 
			filepath.Base(filename)+".html")
		reportFile, err := os.Create(reportPath)
		if err != nil {
			result.Error = fmt.Errorf("error creating report: %w", err)
			return result
		}
		defer reportFile.Close()

		if err := GenerateHTMLReport(stats, reportFile); err != nil {
			result.Error = fmt.Errorf("error generating report: %w", err)
			return result
		}
	}

	// Filter if output directory is specified
	if config.OutputDir != "" {
		input.Seek(0, io.SeekStart) // Reset file pointer
		
		outputPath := filepath.Join(config.OutputDir, 
			"filtered_"+filepath.Base(filename))
		output, err := os.Create(outputPath)
		if err != nil {
			result.Error = fmt.Errorf("error creating output: %w", err)
			return result
		}
		defer output.Close()

		if err := Filter(input, output, config.IsFastq, config.FilterOpts); err != nil {
			result.Error = fmt.Errorf("error filtering: %w", err)
			return result
		}
	}

	return result
}

package skewer

import (
	"fmt"
	"io"
	"sync"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
)

// BatchJob represents a single file processing job.
type BatchJob struct {
	InputFile  string
	OutputFile string
	Index      int
}

// BatchResult contains the result of processing a single job.
type BatchResult struct {
	Job   BatchJob
	Stats *TrimStats
	Error error
}

// ProcessBatch processes multiple FASTQ files in parallel.
func ProcessBatch(jobs []BatchJob, encoding fastq.QualityEncoding, opts TrimOptions, workers int) ([]BatchResult, error) {
	if workers <= 0 {
		workers = 1
	}

	// Create channels for jobs and results
	jobChan := make(chan BatchJob, len(jobs))
	resultChan := make(chan BatchResult, len(jobs))

	// Start worker goroutines
	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for job := range jobChan {
				result := processSingleJob(job, encoding, opts)
				resultChan <- result
			}
		}()
	}

	// Send jobs to workers
	for _, job := range jobs {
		jobChan <- job
	}
	close(jobChan)

	// Wait for all workers to finish
	go func() {
		wg.Wait()
		close(resultChan)
	}()

	// Collect results
	results := make([]BatchResult, 0, len(jobs))
	for result := range resultChan {
		results = append(results, result)
	}

	return results, nil
}

// processSingleJob processes a single FASTQ file.
func processSingleJob(job BatchJob, encoding fastq.QualityEncoding, opts TrimOptions) BatchResult {
	result := BatchResult{Job: job}

	// Open input file
	input, err := iohelper.OpenReader(job.InputFile)
	if err != nil {
		result.Error = fmt.Errorf("error opening input file %s: %w", job.InputFile, err)
		return result
	}
	defer input.Close()

	// Open output file
	output, err := iohelper.OpenWriter(job.OutputFile)
	if err != nil {
		result.Error = fmt.Errorf("error creating output file %s: %w", job.OutputFile, err)
		return result
	}
	defer output.Close()

	// Process file
	stats, err := TrimSingleEnd(input, output, encoding, opts)
	if err != nil {
		result.Error = fmt.Errorf("error processing file %s: %w", job.InputFile, err)
		return result
	}

	result.Stats = stats
	return result
}

// BatchPairedJob represents a paired-end file processing job.
type BatchPairedJob struct {
	InputFile1  string
	InputFile2  string
	OutputFile1 string
	OutputFile2 string
	OutputSingle string
	Index       int
}

// BatchPairedResult contains the result of processing a paired-end job.
type BatchPairedResult struct {
	Job   BatchPairedJob
	Stats *TrimStats
	Error error
}

// ProcessPairedBatch processes multiple paired-end FASTQ files in parallel.
func ProcessPairedBatch(jobs []BatchPairedJob, encoding fastq.QualityEncoding, opts TrimOptions, workers int) ([]BatchPairedResult, error) {
	if workers <= 0 {
		workers = 1
	}

	// Create channels for jobs and results
	jobChan := make(chan BatchPairedJob, len(jobs))
	resultChan := make(chan BatchPairedResult, len(jobs))

	// Start worker goroutines
	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for job := range jobChan {
				result := processPairedJob(job, encoding, opts)
				resultChan <- result
			}
		}()
	}

	// Send jobs to workers
	for _, job := range jobs {
		jobChan <- job
	}
	close(jobChan)

	// Wait for all workers to finish
	go func() {
		wg.Wait()
		close(resultChan)
	}()

	// Collect results
	results := make([]BatchPairedResult, 0, len(jobs))
	for result := range resultChan {
		results = append(results, result)
	}

	return results, nil
}

// processPairedJob processes a single paired-end FASTQ file pair.
func processPairedJob(job BatchPairedJob, encoding fastq.QualityEncoding, opts TrimOptions) BatchPairedResult {
	result := BatchPairedResult{Job: job}

	// Open input files
	input1, err := iohelper.OpenReader(job.InputFile1)
	if err != nil {
		result.Error = fmt.Errorf("error opening input file %s: %w", job.InputFile1, err)
		return result
	}
	defer input1.Close()

	input2, err := iohelper.OpenReader(job.InputFile2)
	if err != nil {
		result.Error = fmt.Errorf("error opening input file %s: %w", job.InputFile2, err)
		return result
	}
	defer input2.Close()

	// Open output files
	output1, err := iohelper.OpenWriter(job.OutputFile1)
	if err != nil {
		result.Error = fmt.Errorf("error creating output file %s: %w", job.OutputFile1, err)
		return result
	}
	defer output1.Close()

	output2, err := iohelper.OpenWriter(job.OutputFile2)
	if err != nil {
		result.Error = fmt.Errorf("error creating output file %s: %w", job.OutputFile2, err)
		return result
	}
	defer output2.Close()

	// Open optional single output file
	var outSingle io.Writer
	if job.OutputSingle != "" {
		f, err := iohelper.OpenWriter(job.OutputSingle)
		if err != nil {
			result.Error = fmt.Errorf("error creating single output file %s: %w", job.OutputSingle, err)
			return result
		}
		defer f.Close()
		outSingle = f
	}

	// Process files
	stats, err := TrimPairedEnd(input1, input2, output1, output2, outSingle, encoding, opts)
	if err != nil {
		result.Error = fmt.Errorf("error processing files %s and %s: %w", job.InputFile1, job.InputFile2, err)
		return result
	}

	result.Stats = stats
	return result
}

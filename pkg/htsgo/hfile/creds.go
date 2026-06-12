package hfile

import (
	"bufio"
	"os"
	"path/filepath"
	"strings"
)

// awsCredentials holds resolved AWS credentials and the region they were found
// with (if any). An empty AccessKey/SecretKey pair means anonymous access.
type awsCredentials struct {
	AccessKey    string
	SecretKey    string
	SessionToken string
	Region       string
}

// resolveAWSCredentials resolves AWS credentials following the precedence used
// by htslib's hfile_s3.c: environment variables first, then the shared
// credentials file. Values present in the environment are not overridden by
// the file. The selected profile comes from AWS_PROFILE or
// AWS_DEFAULT_PROFILE, defaulting to "default".
func resolveAWSCredentials() awsCredentials {
	c := awsCredentials{
		AccessKey:    os.Getenv("AWS_ACCESS_KEY_ID"),
		SecretKey:    os.Getenv("AWS_SECRET_ACCESS_KEY"),
		SessionToken: os.Getenv("AWS_SESSION_TOKEN"),
	}
	if r := os.Getenv("AWS_DEFAULT_REGION"); r != "" {
		c.Region = r
	} else if r := os.Getenv("AWS_REGION"); r != "" {
		c.Region = r
	}

	// If both key parts are already set from the environment we are done; the
	// environment fully wins over the file.
	if c.AccessKey != "" && c.SecretKey != "" {
		return c
	}

	fileCreds, ok := loadSharedCredentials(awsProfile())
	if !ok {
		return c
	}
	if c.AccessKey == "" {
		c.AccessKey = fileCreds.AccessKey
	}
	if c.SecretKey == "" {
		c.SecretKey = fileCreds.SecretKey
	}
	if c.SessionToken == "" {
		c.SessionToken = fileCreds.SessionToken
	}
	if c.Region == "" {
		c.Region = fileCreds.Region
	}
	return c
}

// awsProfile returns the configured AWS profile name, defaulting to "default".
func awsProfile() string {
	if p := os.Getenv("AWS_PROFILE"); p != "" {
		return p
	}
	if p := os.Getenv("AWS_DEFAULT_PROFILE"); p != "" {
		return p
	}
	return "default"
}

// sharedCredentialsPath returns the path of the AWS shared credentials file,
// honouring AWS_SHARED_CREDENTIALS_FILE and otherwise defaulting to
// ~/.aws/credentials.
func sharedCredentialsPath() string {
	if p := os.Getenv("AWS_SHARED_CREDENTIALS_FILE"); p != "" {
		return p
	}
	home, err := os.UserHomeDir()
	if err != nil || home == "" {
		return ""
	}
	return filepath.Join(home, ".aws", "credentials")
}

// loadSharedCredentials reads the requested profile from the shared
// credentials file. It returns the credentials and true on success, or false
// if the file cannot be read or the profile is absent.
func loadSharedCredentials(profile string) (awsCredentials, bool) {
	path := sharedCredentialsPath()
	if path == "" {
		return awsCredentials{}, false
	}
	f, err := os.Open(path)
	if err != nil {
		return awsCredentials{}, false
	}
	defer f.Close()

	var (
		creds   awsCredentials
		inSec   bool
		found   bool
		section string
	)
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") || strings.HasPrefix(line, ";") {
			continue
		}
		if strings.HasPrefix(line, "[") && strings.HasSuffix(line, "]") {
			section = strings.TrimSpace(line[1 : len(line)-1])
			inSec = section == profile
			if inSec {
				found = true
			}
			continue
		}
		if !inSec {
			continue
		}
		key, val, ok := splitINIPair(line)
		if !ok {
			continue
		}
		switch key {
		case "aws_access_key_id":
			creds.AccessKey = val
		case "aws_secret_access_key":
			creds.SecretKey = val
		case "aws_session_token":
			creds.SessionToken = val
		case "region":
			creds.Region = val
		}
	}
	if err := sc.Err(); err != nil {
		return awsCredentials{}, false
	}
	return creds, found
}

// splitINIPair splits a "key = value" INI line into its trimmed, lower-cased
// key and its trimmed value. It reports false when the line has no '='.
func splitINIPair(line string) (key, val string, ok bool) {
	i := strings.IndexByte(line, '=')
	if i < 0 {
		return "", "", false
	}
	key = strings.ToLower(strings.TrimSpace(line[:i]))
	val = strings.TrimSpace(line[i+1:])
	return key, val, true
}

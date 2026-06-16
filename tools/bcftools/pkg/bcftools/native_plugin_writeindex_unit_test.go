// Binary-free unit tests for parseWriteIndexArg (native_plugin_writeindex.go).
// These pin the (fmt, handled, err) mapping for every -W / --write-index
// spelling without touching the upstream binary or building any index.
package bcftools

import "testing"

// TestUnitWriteIndexArg covers the bare flag (CSI default), explicit csi/tbi
// suffixes (case-insensitive), the empty suffix (CSI), a non-write-index token
// (handled=false), and an unknown suffix (error).
func TestUnitWriteIndexArg(t *testing.T) {
	tests := []struct {
		name        string
		arg         string
		wantFmt     writeIndexFmt
		wantHandled bool
		wantErr     bool
	}{
		{"bare short flag selects CSI", "-W", writeIndexCSI, true, false},
		{"bare long flag selects CSI", "--write-index", writeIndexCSI, true, false},
		{"short csi", "-W=csi", writeIndexCSI, true, false},
		{"short tbi", "-W=tbi", writeIndexTBI, true, false},
		{"long csi", "--write-index=csi", writeIndexCSI, true, false},
		{"long tbi", "--write-index=tbi", writeIndexTBI, true, false},
		{"uppercase suffix is case-insensitive", "-W=TBI", writeIndexTBI, true, false},
		{"mixed-case suffix", "--write-index=Csi", writeIndexCSI, true, false},
		{"empty short suffix is CSI default", "-W=", writeIndexCSI, true, false},
		{"empty long suffix is CSI default", "--write-index=", writeIndexCSI, true, false},
		{"unrelated token is not handled", "-O", writeIndexOff, false, false},
		{"unknown suffix errors", "-W=bogus", writeIndexOff, true, true},
		{"unknown long suffix errors", "--write-index=zzz", writeIndexOff, true, true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotFmt, gotHandled, err := parseWriteIndexArg(tt.arg)
			if (err != nil) != tt.wantErr {
				t.Fatalf("parseWriteIndexArg(%q) err = %v, wantErr %v", tt.arg, err, tt.wantErr)
			}
			if gotFmt != tt.wantFmt {
				t.Errorf("parseWriteIndexArg(%q) fmt = %v, want %v", tt.arg, gotFmt, tt.wantFmt)
			}
			if gotHandled != tt.wantHandled {
				t.Errorf("parseWriteIndexArg(%q) handled = %v, want %v", tt.arg, gotHandled, tt.wantHandled)
			}
		})
	}
}

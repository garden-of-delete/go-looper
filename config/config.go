package config

// Config holds the global configuration values for the application
type Config struct {
	NucleationFreeEnergy *float64
	SuperhelicityDomain  *int
	AutoDomainSize       bool
	SuperhelicalDensity  *float64
	MinRLoopLength       *int
	Reverse              *bool
	Complement           *bool
	Unconstrained        *bool
	Invert               bool
	Dump                 bool
	Circular             bool
	Residuals            bool
	LocalAverageEnergy   bool
	Homopolymer          *float64
	InfileName           string
	OutfileName          string
}

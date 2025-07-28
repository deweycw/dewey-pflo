# PFLOTRAN Performance Monitor

A comprehensive performance monitoring tool designed to identify bottlenecks in PFLOTRAN reactive transport simulations and guide hardware optimization decisions.

## What This Script Does

The PFLOTRAN Performance Monitor provides detailed analysis of:
- **Memory bandwidth utilization** and cache performance
- **CPU efficiency** across different core types (P-cores vs E-cores on Intel 13th gen)
- **Storage I/O patterns** and bottlenecks
- **System resource usage** over time
- **Specific recommendations** for hardware upgrades

**Primary Goal**: Determine whether your PFLOTRAN performance is limited by memory bandwidth, storage speed, or CPU capacity.

## Prerequisites

### System Requirements
- Linux system with `perf` tools installed
- Intel 13th Gen CPU (script optimized for hybrid P+E core architecture)
- `kernel.perf_event_paranoid = 1` (or lower) for performance monitoring

### Required Software
```bash
# Ubuntu/Debian
sudo apt install linux-tools-common linux-tools-generic python3

# RHEL/CentOS/Rocky
sudo yum install perf python3

# Check installation
perf --version
python3 --version
```

### System Configuration
```bash
# Enable performance monitoring (required once)
echo 'kernel.perf_event_paranoid = 1' | sudo tee -a /etc/sysctl.conf
sudo sysctl -p

# Verify setting
cat /proc/sys/kernel/perf_event_paranoid
# Should show: 1
```

## Installation and Setup

1. **Download the script**:
   ```bash
   wget https://your-script-location/pflotran_working_monitor.sh
   # or copy the script to your system
   ```

2. **Make it executable**:
   ```bash
   chmod +x pflotran_working_monitor.sh
   ```

3. **Test perf functionality** (optional but recommended):
   ```bash
   wget https://your-script-location/perf_test_pflotran.sh
   chmod +x perf_test_pflotran.sh
   ./perf_test_pflotran.sh
   ```

## How to Use

### Basic Usage
```bash
./pflotran_working_monitor.sh mpirun -n 16 ./pflotran -input_prefix copper_leach
```

### With NUMA Optimization
```bash
./pflotran_working_monitor.sh numactl --interleave=all mpirun -n 16 --bind-to core ./pflotran -input_prefix copper_leach
```

### Examples
```bash
# Small test run
./pflotran_working_monitor.sh mpirun -n 4 ./pflotran -input_prefix test_case

# Production run with 32 processes
./pflotran_working_monitor.sh mpirun -n 32 ./pflotran -input_prefix large_simulation

# With specific MPI options
./pflotran_working_monitor.sh mpirun -n 16 --bind-to core --map-by core ./pflotran -input_prefix mycase
```

## Output Files and Directory Structure

The script creates a timestamped directory with comprehensive monitoring data:

```
pflotran_perf_20250727_143245/
├── performance_analysis.txt    # 📋 START HERE - Main analysis and recommendations
├── final_report.txt           # 📊 Executive summary
├── perf_raw.log              # 🔧 Raw performance counter data
├── memory.log                # 💾 Memory usage over time
├── iostat.log                # 💽 Storage I/O statistics  
├── vmstat.log                # ⚙️ CPU and system resource usage
├── memory_summary.txt        # 📈 Memory usage summary
├── storage_summary.txt       # �� Storage performance summary
└── analyze_intel_hybrid.py   # 🐍 Analysis script (auto-generated)
```

### Key Files to Review

| File | Purpose | When to Check |
|------|---------|---------------|
| `performance_analysis.txt` | **Main results and recommendations** | Always - start here |
| `final_report.txt` | Quick overview and file guide | First-time users |
| `memory_summary.txt` | Memory usage patterns | If memory issues suspected |
| `storage_summary.txt` | I/O performance summary | If storage bottlenecks suspected |

## Understanding the Output

### 1. Performance Analysis Structure

The main analysis (`performance_analysis.txt`) contains:

#### **Core Utilization**
```
Core Type Utilization:
  P-cores (Performance): ✓ Active
  E-cores (Efficiency):  ✓ Active
```
- Shows whether both P-cores and E-cores are being used
- Both should be active for optimal performance

#### **Cache Performance**
```
Cache Performance (Combined P+E cores):
  Cache misses: 96,822,254,944
  Cache references: 114,179,070,991
  Cache miss ratio: 84.80%
  🔴 CRITICAL - Extreme cache miss rate indicates severe memory bandwidth bottleneck
```

#### **CPU Efficiency**
```
CPU Efficiency (Combined P+E cores):
  Instructions per cycle: 1.95
  ✓ GOOD - Decent CPU utilization
```

#### **Per-Core Analysis**
```
Per-Core Type Analysis:
  P-cores IPC: 2.72
  E-cores IPC: 0.95
  Workload distribution: 78.9% P-cores, 21.1% E-cores
```

### 2. Performance Assessment Categories

#### 🟢 **GOOD PERFORMANCE**
- Cache miss ratio: <10%
- Instructions per cycle: >1.5
- **Action**: System is well-balanced, no urgent upgrades needed

#### 🟡 **MODERATE ISSUES**
- Cache miss ratio: 10-25%
- Instructions per cycle: 1.0-1.5
- **Action**: Monitor performance, consider upgrades if budget allows

#### 🔴 **SEVERE BOTTLENECK**
- Cache miss ratio: >50%
- Instructions per cycle: <1.0
- **Action**: Immediate hardware upgrade recommended

#### 🔴 **CRITICAL BOTTLENECK**
- Cache miss ratio: >80%
- **Action**: System performance severely degraded, urgent upgrade required

### 3. Cache Miss Ratio Interpretation

| Cache Miss Ratio | Performance Level | Meaning | Action Required |
|------------------|-------------------|---------|-----------------|
| <5% | ✅ Excellent | Data fits in cache, optimal performance | None |
| 5-10% | ✅ Good | Acceptable cache efficiency | None |
| 10-15% | ⚠️ Fair | Some cache pressure | Monitor |
| 15-25% | ⚠️ Poor | Notable cache inefficiency | Consider memory upgrade |
| 25-50% | 🔴 Severe | Major memory bandwidth issues | Upgrade recommended |
| >50% | 🔴 Critical | Extreme memory bandwidth bottleneck | **Urgent upgrade required** |

### 4. Specific Recommendations

The script provides tailored recommendations based on detected bottlenecks:

#### **Memory Bandwidth Bottleneck**
```
🔴 SEVERE MEMORY BANDWIDTH BOTTLENECK
   Immediate upgrade to DDR5-5600+ recommended for substantial performance gains
   Expected improvement: 60-80% faster simulation times
```
**Action**: Upgrade from DDR4 to DDR5 memory

#### **Storage Bottleneck**
```
🟡 STORAGE I/O LIMITATIONS
   Consider PCIe 4.0 NVMe SSD for improved I/O performance
   Expected improvement: 20-30% faster startup and output writing
```
**Action**: Upgrade to faster storage

#### **CPU Limitations**
```
⚠️ CPU UTILIZATION ISSUES
   Consider higher core count processor for better parallel performance
```
**Action**: Upgrade CPU or optimize MPI configuration

### 5. Memory Bandwidth Impact Analysis

For severe cases (>50% cache miss rate), the script calculates specific performance impact:

```
Memory Bandwidth Impact Analysis:
  With 84.8% cache miss rate:
  Current DDR4-3200: ~255 cycles average per memory access
  With DDR5-5600: ~153 cycles average per memory access
  Expected improvement: 1.7x faster with DDR5 upgrade
```

This quantifies the exact performance benefit of a memory upgrade.

## Troubleshooting

### Common Issues

#### 1. **Permission Denied for perf**
```bash
# Error: Access to performance monitoring and observability operations is limited
sudo sysctl kernel.perf_event_paranoid=1
echo 'kernel.perf_event_paranoid = 1' | sudo tee -a /etc/sysctl.conf
```

#### 2. **perf Events Not Available**
```bash
# Error: event syntax error: 'cache-misses'
# Run diagnostic test
./perf_test_pflotran.sh
# Use only working events from diagnostic output
```

#### 3. **MPI + perf Compatibility Issues**
```bash
# Try alternative MPI + perf ordering
mpirun -n 16 perf stat -e cache-misses ./pflotran -input_prefix mycase
# Instead of: perf stat -e cache-misses mpirun -n 16 ./pflotran
```

#### 4. **Empty Output Files**
```bash
# Check if background processes started properly
ps aux | grep iostat
ps aux | grep vmstat
# Restart with verbose output to debug
```

### Performance Monitoring Not Working?

If perf monitoring fails, the script falls back to basic system monitoring:
- Timing with `/usr/bin/time -v`
- Memory usage tracking
- I/O statistics with `iostat`
- CPU usage with `vmstat`

This still provides valuable bottleneck identification information.

## Optimization Workflow

### Phase 1: Identify Bottleneck
1. Run monitoring script with typical PFLOTRAN workload
2. Review `performance_analysis.txt`
3. Identify primary bottleneck (memory/storage/CPU)

### Phase 2: Quick Optimizations
```bash
# Try NUMA optimization
numactl --interleave=all mpirun -n 16 ./pflotran -input_prefix mycase

# Test matrix reordering in PFLOTRAN input file
MATRIX_ORDERING RCM

# Run monitoring again to measure improvement
```

### Phase 3: Hardware Upgrade
Based on analysis results:
- **Memory bottleneck**: Upgrade DDR4 → DDR5
- **Storage bottleneck**: Upgrade to PCIe 4.0 NVMe SSD  
- **CPU bottleneck**: Higher core count or faster CPU

### Phase 4: Validate Improvements
Re-run monitoring script after changes to quantify performance gains.

## Example Interpretation

### Scenario 1: Memory Bandwidth Limited
```
Cache miss ratio: 84.80%
�� SEVERE MEMORY BANDWIDTH BOTTLENECK
Expected improvement: 1.7x faster with DDR5 upgrade
```
**Interpretation**: DDR4 memory cannot keep up with PFLOTRAN's data demands. Upgrade to DDR5 for major performance improvement.

### Scenario 2: Storage Limited
```
Cache miss ratio: 8.2%
Storage utilization: 95%
🟡 STORAGE I/O LIMITATIONS  
```
**Interpretation**: Memory performance is good, but storage is saturated. Upgrade to faster SSD.

### Scenario 3: Well Balanced System
```
Cache miss ratio: 6.1%
Instructions per cycle: 2.3
🟢 GOOD PERFORMANCE
```
**Interpretation**: System is well-optimized for current workload. No immediate upgrades needed.

## Quick Reference Commands

```bash
# View main analysis
cat pflotran_perf_*/performance_analysis.txt

# Check memory usage summary  
cat pflotran_perf_*/memory_summary.txt

# Check storage performance
cat pflotran_perf_*/storage_summary.txt

# View all results
ls -la pflotran_perf_*/

# Find peak memory usage
grep -o '[0-9.]*Gi' pflotran_perf_*/memory.log | sed 's/Gi//' | sort -n | tail -1

# Check storage bottlenecks  
awk '{if($7>80) print "High utilization: " $0}' pflotran_perf_*/iostat.log
```

## Support and Contributing

### Getting Help
- Review this README and troubleshooting section
- Run the diagnostic script: `./perf_test_pflotran.sh`  
- Check system requirements and prerequisites

### Improving the Script
The monitoring script is designed to be extensible. Key areas for enhancement:
- Additional CPU architectures (AMD, older Intel)
- More sophisticated PFLOTRAN-specific analysis
- Integration with PFLOTRAN's built-in profiling
- Automated optimization recommendations

---

**Remember**: The goal is to identify whether memory bandwidth, storage speed, or CPU capacity is limiting your PFLOTRAN performance, then make targeted hardware upgrades for maximum performance improvement per dollar spent.

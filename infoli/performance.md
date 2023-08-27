Time in microseconds for different network sizes (horizontal) and code versions
(vertical):

|         | 10      | 100       | 250       | 500       | 750        | 1000       |
|---------|---------|-----------|-----------|-----------|------------|------------|
| serial  | 73,022  | 785,708   | 3,022,811 | 9,787,185 | 23,091,656 | 38,038,550 |
|---------|---------|-----------|-----------|-----------|------------|------------|
| 31796a6 | 850,326 | 1,083,953 | 1,236,045 | 1,687,585 | 1,957,936  | 2,316,141  |
| SPEEDUP | 0.1X    | 0.7X      | 2.4X      | 5.8X      | 11.8X      | 16.4X      |
|---------|---------|-----------|-----------|-----------|------------|------------|
| 8af2595 | 821,856 | 1,200,881 | 1,262,908 | 1,429,628 | 1,696,038  | 2,058,642  |
| SPEEDUP | 0.1X    | 0.7X      | 2.4X      | 6.8X      | 13.6X      | 18.5X      |

The serial version contains logic that is O(N^2) in the network size, which can be
clearly seen here from the deterioration in runtime as the network size grows. The
parallel version(s) on the other hand, effectively parallelizes the O(N^2) workload such
that it becomes increasingly faster than the serial version as the network size
increases. For small network sizes, the serial version outperforms the parallel version,
most likely due to additional overhead for copying to and from GPU memory.

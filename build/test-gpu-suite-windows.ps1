# Runs a GPU program on every work unit of dataFiles against its truth file (Windows).
#   -exe       the program: a file in <work>\bin, or a full path       -tag   a name for the run
#   -df        upstream's Test\dataFiles                               -work  where the run directories go
#   -initData  an init_data.xml for each run (default: <work>\gpu\init_data.xml if it exists)
#   -envs      variables set during the runs, "NAME=value;..." (e.g. NF_OCL_ANY=1 or NF_CUDA_ANY=1 without an init_data.xml)
# The defaults are the original author's machine.
param([string]$exe, [string]$tag,
      [string]$df = 'C:\Users\Alp\AppData\Local\Temp\claude\C--Users-Alp-Documents-Default-Project\ac77d448-3f34-4416-ac54-fafab05b597e\scratchpad\src\Test\dataFiles',
      [string]$work = 'C:\Users\Alp\AppData\Local\Temp\nfw', [string]$initData = '', [string]$envs = '')
$b = $work
if (-not $initData -and (Test-Path "$b\gpu\init_data.xml")) { $initData = "$b\gpu\init_data.xml" }
$prog = if (Test-Path $exe) { (Resolve-Path $exe).Path } else { "$b\bin\$exe" }
foreach ($kv in ($envs -split ';' | Where-Object { $_ })) { $k,$v = $kv -split '=',2; Set-Item -Path "Env:$k" -Value $v }
$ok=0; $bad=0; $tot=0.0
Get-ChildItem "$df\*DS*_Grp*.dat" | ForEach-Object {
  $wu = $_.BaseName
  $d="$b\suite2_$tag`_$wu"; Remove-Item -Recurse -Force $d -ErrorAction SilentlyContinue; New-Item -ItemType Directory $d | Out-Null
  if ($initData) { Copy-Item $initData "$d\init_data.xml" }; Copy-Item "$df\$wu.dat" "$d\in"
  $t0 = Get-Date
  $p = Start-Process -FilePath $prog -WorkingDirectory $d -PassThru -Wait -NoNewWindow
  $s = ((Get-Date) - $t0).TotalSeconds; $tot += $s
  $o = Get-Content "$d\out" -ErrorAction SilentlyContinue | Where-Object { $_ -notmatch 'Elapsed Time' }
  $tr = Get-Content "$df\wu_$wu.out.truth" | Where-Object { $_ -notmatch 'Elapsed Time' }
  $r = if ($o -and ((Compare-Object $o $tr) -eq $null) -and $p.ExitCode -eq 0) { 'MATCH' } else { 'DIFF' }
  if ($r -eq 'MATCH') { $ok++ } else { $bad++ }
  "$tag $wu $r exit=$($p.ExitCode) wall=$([math]::Round($s,1)) s"
}
"SUITE $tag : MATCH $ok  DIFF $bad  total $([math]::Round($tot,0)) s"
foreach ($kv in ($envs -split ';' | Where-Object { $_ })) { Remove-Item "Env:$(($kv -split '=',2)[0])" }

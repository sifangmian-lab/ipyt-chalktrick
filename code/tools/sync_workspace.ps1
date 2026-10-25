[CmdletBinding()]
param(
    [string]$LocalPath = 'D:\vibelab\iypt',
    [string]$Branch = 'develop'
)

$ErrorActionPreference = 'Stop'
$Repository = 'https://github.com/sifangmian-lab/ipyt-chalktrick.git'

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw '未检测到 Git。请先安装 Git for Windows，再重新运行本脚本。'
}

$FullPath = [System.IO.Path]::GetFullPath($LocalPath)
$GitDir = Join-Path $FullPath '.git'
$SnapshotCommit = $null

if (-not (Test-Path -LiteralPath $FullPath)) {
    $Parent = Split-Path -Parent $FullPath
    New-Item -ItemType Directory -Force -Path $Parent | Out-Null

    & git clone --branch $Branch --single-branch $Repository $FullPath
    if ($LASTEXITCODE -ne 0) {
        throw "首次克隆失败；目标路径：$FullPath"
    }
}
elseif (-not (Test-Path -LiteralPath $GitDir)) {
    $ExistingItems = @(Get-ChildItem -LiteralPath $FullPath -Force)

    if ($ExistingItems.Count -eq 0) {
        & git clone --branch $Branch --single-branch $Repository $FullPath
        if ($LASTEXITCODE -ne 0) {
            throw "首次克隆失败；目标路径：$FullPath"
        }
    }
    else {
        Write-Host "检测到已有文件，将先建立本地快照，再并入 origin/$Branch。"

        $ConfiguredName = (& git -C $FullPath config user.name 2>$null)
        $ConfiguredEmail = (& git -C $FullPath config user.email 2>$null)
        if ([string]::IsNullOrWhiteSpace($ConfiguredName) -or [string]::IsNullOrWhiteSpace($ConfiguredEmail)) {
            throw '建立本地快照需要 Git 用户名和邮箱。请先配置 user.name 与 user.email；脚本尚未修改该目录。'
        }

        & git ls-remote $Repository "refs/heads/$Branch" | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "无法访问远端 $Branch 分支；脚本尚未修改该目录。"
        }

        & git -C $FullPath init -b $Branch
        if ($LASTEXITCODE -ne 0) {
            throw '无法在现有目录初始化 Git。'
        }

        & git -C $FullPath remote add origin $Repository
        if ($LASTEXITCODE -ne 0) {
            throw '无法登记 origin 远程仓库。'
        }

        & git -C $FullPath fetch origin $Branch
        if ($LASTEXITCODE -ne 0) {
            throw '获取远端更新失败。'
        }

        $IgnoredItems = @(& git -C $FullPath ls-files --others --ignored --exclude-standard)
        if ($IgnoredItems.Count -gt 0) {
            Write-Warning "有 $($IgnoredItems.Count) 个被忽略路径不会进入 Git 快照；重要文件请先解除忽略或另行备份。"
        }

        & git -C $FullPath add --all
        if ($LASTEXITCODE -ne 0) {
            throw '无法暂存现有文件。'
        }

        & git -C $FullPath commit -m '合并远端前保存本地现有文件'
        if ($LASTEXITCODE -ne 0) {
            throw '无法建立本地快照。请检查文件是否均被忽略，以及 Git 用户名和邮箱配置。'
        }

        $SnapshotCommit = (& git -C $FullPath rev-parse HEAD).Trim()

        & git -C $FullPath merge --allow-unrelated-histories --no-edit "origin/$Branch"
        if ($LASTEXITCODE -ne 0) {
            if (Test-Path -LiteralPath (Join-Path $GitDir 'MERGE_HEAD')) {
                & git -C $FullPath merge --abort | Out-Null
            }
            throw "自动合并未完成，可能同名文件内容不同。未强制覆盖；本地快照为 $SnapshotCommit。处理冲突后可重新合并 origin/$Branch。"
        }

        & git -C $FullPath branch "--set-upstream-to=origin/$Branch" $Branch
        if ($LASTEXITCODE -ne 0) {
            throw '文件已合并，但设置上游分支失败；请检查 git status。'
        }
    }
}
else {
    $RemoteUrl = (& git -C $FullPath remote get-url origin 2>$null)
    if ($LASTEXITCODE -ne 0) {
        throw '现有工作区缺少 origin 远程仓库，脚本未做修改。'
    }
    if ($RemoteUrl -notmatch 'github\.com[:/]sifangmian-lab/ipyt-chalktrick(?:\.git)?$') {
        throw "现有工作区的 origin 不是本项目，脚本已停止：$RemoteUrl"
    }

    & git -C $FullPath fetch origin $Branch
    if ($LASTEXITCODE -ne 0) {
        throw '获取远端更新失败。'
    }

    & git -C $FullPath checkout $Branch
    if ($LASTEXITCODE -ne 0) {
        throw "无法切换到 $Branch 分支；请先处理本地未提交修改。"
    }

    & git -C $FullPath pull --ff-only origin $Branch
    if ($LASTEXITCODE -ne 0) {
        throw '无法快进同步。可能存在本地提交或分支分歧；脚本未强制覆盖，请检查 git status。'
    }
}

Write-Host ''
Write-Host "工作区路径：$FullPath"
Write-Host "跟踪分支：$Branch"
if ($SnapshotCommit) {
    Write-Host "合并前本地快照：$SnapshotCommit"
}
& git -C $FullPath status --short --branch
Write-Host ''
Write-Host '同步完成。脚本不会自动提交或推送本地修改。'

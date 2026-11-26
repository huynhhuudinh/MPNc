# MPNc model: tái lập Application và mô phỏng soft-closure

Kho mã này tái lập phần **Application** (Mục 7) của bài báo *Identifiability
and Observational Equivalence in Soft-Closure Multinomial N-Mixture Models* trên
dữ liệu Chestnut-sided Warbler (`cswa`), và cung cấp CLI mô phỏng (Mục 6) cho
bốn mô hình MPN0, MPN1, MPNc và CRTE.

Toàn bộ phân tích được viết bằng Python; không cần cài R. Tệp nguồn `.rda` được
đọc bằng gói `rdata` nếu đã cài, nếu không thì bằng bộ giải mã không phụ thuộc
`rda_reader.py`; kiểm thử xác nhận hai cách cho cùng dữ liệu. Báo cáo chi tiết
nằm trong [`APPLICATION_REPORT.md`](APPLICATION_REPORT.md).

## Trạng thái kết quả đi kèm

Kết quả gốc được tạo bằng Python 3.12.14 (NumPy 2.3.5, SciPy 1.17.0) và đã
được chạy lại bằng Python 3.13.16 (NumPy 2.5.3, SciPy 1.18.1). Kết quả trong
`results/` của gói này là từ lần chạy lại gần nhất, bằng Python 3.12.3
(NumPy 2.2.6, SciPy 1.16.3) với bộ giải mã `rda_reader`
(`results/run_manifest.json`), và trùng với các lần chạy trước đến ít nhất 6
chữ số thập phân. Cả 8 phép khớp đều thành công; **77/77** con số của
Mục 7 trong bản thảo (ước lượng, logLik, AIC, trọng số, khoảng tin cậy, chẩn
đoán bậc ba, kiểm định composition, tổng diện tích) khớp sau làm tròn
(`results/paper_comparison.csv`).

| Phân tích | Mô hình tốt nhất theo AIC | AIC | Ước lượng abundance | `r` mỗi phút | Tham số phụ thuộc |
| --- | --- | ---: | ---: | ---: | ---: |
| Offset diện tích | MPNc | 353.666008 | 1.083596 cá thể/ha | 0.373167 | `c = 0.456632` |
| Chỉ hệ số chặn | MPNc | 356.120210 | 0.589969 cá thể/điểm/lần | 0.371773 | `c = 0.456554` |

Đây là kết quả trên dữ liệu thực, không phải smoke test mô phỏng. Xem
`results/verification_summary.json` để kiểm tra hội tụ và
`results/paper_comparison.csv` để xem từng phép đối chiếu.

## Cài đặt

Các phiên bản thư viện dùng cho lần chạy đã giao được cố định trong
`requirements.txt`. `rdata` và `xarray` là tùy chọn: thiếu chúng, mã tự dùng
`rda_reader.py`.

```bash
cd MPNc_model          # thư mục giải nén của gói
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Trên Windows, dùng `.venv\Scripts\activate` để kích hoạt môi trường. Các script
xác định đường dẫn mặc định theo vị trí của chính script nên cũng có thể được
gọi từ thư mục khác bằng đường dẫn tuyệt đối.

## Chạy lại Application

Workflow đầy đủ gồm kiểm thử, chuẩn bị dữ liệu, khớp mô hình và dựng báo cáo:

```bash
python -m unittest -v test_application.py
python run_application.py --starts 8
python build_report.py
```

`run_application.py` chạy hai đặc tả (`area` và `intercept`) cho cả bốn mô hình.
Các tùy chọn hiện có:

```text
--data-dir PATH   thư mục chứa dữ liệu; mặc định là data/
--out PATH        thư mục kết quả; mặc định là results/
--download        tải lại ba tệp nguồn đã ghim và kiểm tra SHA-256
--prepare-only    chỉ kiểm tra và xuất dữ liệu chuẩn hóa, không khớp mô hình
--starts N        số điểm khởi đầu; mặc định là 8 và phải dương
```

Ví dụ giữ nguyên bộ kết quả đã giao:

```bash
python run_application.py --out results_new --starts 8
python build_report.py --results results_new
```

Lần chạy mới chỉ ghi đè các tệp đầu ra trùng tên trong thư mục đã chọn; không
xóa các tệp khác ở đó. `build_report.py --results ...` đọc kết quả từ thư mục
được chỉ định, tạo các bảng kiểm tra trong chính thư mục đó và luôn cập nhật
`APPLICATION_REPORT.md` tại thư mục dự án.

### Dữ liệu

`data/raw/` chứa `cswa.rda`, tài liệu `data_cswa.Rd`, metadata gói
`DESCRIPTION` và `provenance.json`. Nguồn được ghim tại một commit của bản sao
CRAN `AHMbook`; hash SHA-256 được kiểm tra trước mỗi lần đọc. Không cần mạng khi
dùng các tệp đã có.

```bash
# Tải lại nguồn đã ghim, kiểm tra hash rồi chạy phân tích
python run_application.py --download --starts 8

# Chỉ xác minh nguồn và tái tạo data/*.csv, data/*.npz, data_audit.json
python run_application.py --prepare-only
```

Dữ liệu có 43 điểm, 3 lần khảo sát và các khoảng removal dài 2, 3, 5 phút. Có
128 lượt khảo sát quan sát được; toàn bộ lần 1 của điểm `darty` bị thiếu và được
loại khỏi likelihood. Trong `counts_long.csv`, hàng của lần thiếu vẫn được giữ
với `count` trống và `observed=False`; trong `analysis_data.npz`, giá trị này vẫn
là `NaN`, không bị đổi thành 0. Covariate diện tích được ghép theo `point_id`,
không theo vị trí hàng.

Các tệp dữ liệu đã chuẩn hóa:

| Tệp | Nội dung |
| --- | --- |
| `data/counts_long.csv` | Dữ liệu đếm dạng long, gồm cờ quan sát và diện tích điểm |
| `data/site_covariates.csv` | Covariate theo điểm, đã căn chỉnh bằng định danh |
| `data/spot_mapping.csv` | Dữ liệu spot-mapping theo patch |
| `data/analysis_data.npz` | Mảng `y`, `area_ha` và `point_id` dùng trong phân tích |
| `data/data_audit.json` | Trục dữ liệu, missingness, phép ghép, đơn vị và các tổng kiểm tra |

### Mô hình và đơn vị

- **MPN0**: abundance ở các lần khảo sát độc lập; đây là working likelihood khi
  dữ liệu thực tế có phụ thuộc. Sai số chuẩn chính dùng sandwich gom cụm theo
  điểm.
- **MPN1**: abundance đóng hoàn toàn qua các lần khảo sát.
- **MPNc**: thành phần Poisson chung có mean `c * mu_i`, thành phần riêng có mean
  `(1 - c) * mu_i`.
- **CRTE**: superpopulation Poisson có mean `mu_i / phi`; availability độc lập
  qua lần khảo sát với xác suất `phi`. `mu_i` là mean hiện diện trong một lần.

Trong phân tích `area`, `mu_i = rho * plotArea_i` và `rho` có đơn vị cá thể/ha.
Trong phân tích `intercept`, `mu_i = mu` và đơn vị là cá thể/điểm/lần khảo sát.
`r` luôn là xác suất phát hiện mỗi phút. Cột `dep` là `c` với MPNc, `phi` với
CRTE, và được cố định lần lượt ở 0/1 cho MPN0/MPN1.

Likelihood được tính trên log-scale, giữ các hằng số factorial/multinomial để
tái lập logLik và AIC tuyệt đối. Tổng trên abundance ẩn dùng ngưỡng đuôi và một
phép kiểm tra chặt hơn; các giới hạn số nằm trong `FitConfig` của
`softclosure.py`.

### Đầu ra Application

| Tệp | Nội dung |
| --- | --- |
| `results/model_results.csv` | Ước lượng, logLik, AIC, trọng số, SE, CI và chẩn đoán của 8 phép khớp |
| `results/area_*.json` | Chi tiết từng mô hình cho đặc tả offset diện tích |
| `results/intercept_*.json` | Chi tiết từng mô hình cho đặc tả chỉ hệ số chặn |
| `results/site_loglik.csv` | Đóng góp log-likelihood theo điểm và mô hình |
| `results/third_order.json` | Chẩn đoán residual bậc ba theo diện tích và bản pooled để audit |
| `results/composition_gof.json` | Kiểm định độ phù hợp của composition removal (X², G², p bootstrap) |
| `results/run_manifest.json` | Phiên bản phần mềm, cấu hình, nguồn dữ liệu và hash mã đã chạy |
| `results/paper_comparison.csv` | Đối chiếu 77 con số của Mục 7 với dung sai làm tròn |
| `results/confidence_intervals.csv` | CI log/logit và Wald trên thang tự nhiên, kèm loại covariance |
| `results/verification_summary.json` | Tổng hợp hội tụ, độ ổn định số và phép đối chiếu |
| `results/tests.log`, `results/run.log` | Nhật ký kiểm thử và lần chạy đã giao |

Trong JSON/CSV, tiền tố `mu_` từ mô-đun số dùng chung chỉ abundance đầu tiên:
cần đọc là `rho` với phân tích `area` và `mu` với `intercept`.
`confidence_intervals.csv` đã đổi tên tham số theo estimand. Trường `primary_se`
cho biết covariance được chọn làm kết quả chính.

## Chạy mô phỏng

`softclosure.py` vừa cung cấp các hàm số được `run_application.py` tái sử dụng,
vừa là CLI mô phỏng độc lập cho thiết kế có các khoảng removal bằng nhau. Các
suite hiện có là `smoke` (9 kịch bản), `core` (24), `reverse` (24), `small` (6),
`stress` (6) và `full` (60).

Phần tổng hợp mô phỏng tạo biểu đồ PNG/PDF và cần thêm `matplotlib` (không nằm
trong `requirements.txt` của workflow Application):

```bash
python -m pip install matplotlib
python softclosure.py run --suite smoke --reps 10 --workers 2 --out simulation_smoke
```

Nếu bỏ `--reps`, mặc định là 10 cho `smoke` và 500 cho các suite khác. Có thể
dùng `--seed`, `--starts`, `--maxiter` và `--no-intervals`; `--workers` phải từ
1 đến 8. Mỗi replicate được ghi riêng dưới `OUT/replicates/`, cho phép tiếp tục
một lần chạy dở nếu fingerprint của code, seed, scenario và cấu hình vẫn khớp.
Không dùng lại cùng thư mục output khi thay một trong các thành phần đó.

Tạo lại bảng tổng hợp và hình từ các replicate đã hoàn tất, rồi dựng các bảng và
hình dùng trong bài (Bảng tỷ lệ khớp, bảng chọn mô hình, ô sai đặc tả, hai hình
`full_misspecification` và `full_model_selection`) vào `OUT/paper/`:

```bash
python softclosure.py report --out simulation_smoke
python paper_simulation_outputs.py --out simulation_smoke
```

Lần chạy đầy đủ của bài: `python softclosure.py run --suite full --out simulation_full`
(seed mặc định 20260910, 500 replicate, 3 điểm khởi đầu). Khoảng 1,5–4 giây cho
mỗi bộ dữ liệu trên một lõi, tức vài giờ với 8 worker.

Đầu ra mô phỏng gồm `manifest.json`, `fits.csv`, `parameter_summary.csv`,
`model_selection.csv`, `selection_summary.csv`, `third_order.csv`,
`third_order_summary.csv`, `run_summary.json` và các hình `rmse_*` ở định dạng
PNG/PDF. Smoke suite chỉ dùng để kiểm tra đường chạy, không phải bằng chứng mô
phỏng đầy đủ.

## Cấu trúc mã

| Tệp | Vai trò |
| --- | --- |
| `run_application.py` | Xác minh nguồn, chuẩn hóa dữ liệu, khớp 8 mô hình, chẩn đoán bậc ba, kiểm định composition |
| `rda_reader.py` | Bộ giải mã `.rda` (RDX2/RDX3, XDR) không phụ thuộc gói ngoài |
| `softclosure.py` | Likelihood/lõi số dùng chung và CLI mô phỏng có thể resume |
| `paper_simulation_outputs.py` | Dựng bảng và hình mô phỏng của bài từ thư mục kết quả |
| `test_application.py` | 10 kiểm thử: xác suất removal, trục dữ liệu, missingness, bốn likelihood so với tổng multinomial gốc, bộ giải mã độc lập, kiểm định composition, profile working |
| `build_report.py` | Đối chiếu 77 con số với bài, dựng bảng CI, verification summary và báo cáo Markdown |
| `PACKAGE_MANIFEST.json` | Snapshot kích thước và SHA-256 của các tệp được giao, gồm cả `simulation_full/` (trừ 30.000 tệp trong `simulation_full/replicates/`, chỉ ghi số lượng và tổng dung lượng) |
| `LICENSE` | Giấy phép GPL-3.0-or-later |

## Phạm vi và diễn giải

Workflow Application chỉ tái lập hai đặc tả intercept/area; chưa thêm covariate
thời gian, ngày, người quan sát, habitat, hiệu ứng patch hay overdispersion. Dữ
liệu không có định danh cá thể nên `c` không trực tiếp chứng minh cùng một cá
thể quay lại. Phân tích xem 43 điểm là độc lập như bài báo, dù các điểm thuộc 38
patch; giả định này chưa được kiểm tra ở đây. Trọng số AIC không phải xác suất
một cơ chế sinh học là đúng.

Nguồn dữ liệu ghi công David King và cộng sự; tài liệu nguồn dẫn Chandler, King
& Chandler (2009), Chandler, Royle & King (2011), và Kéry & Royle (2021). Gói
AHMbook khai báo GPL (>= 3); khi phân phối lại cần giữ các tệp và provenance
trong `data/raw/`. Mã Python của dự án được cung cấp theo GPL-3.0-or-later; xem
`LICENSE`.

# Kết quả tái lập phần Application

Đã chạy **8 phép khớp: 4 mô hình × 2 cách mô hình hóa abundance** trên dữ liệu thực. 8/8 phép khớp thành công. Đã đối chiếu 77 con số in trong Mục 7 của bản thảo (ước lượng, logLik, AIC, trọng số, khoảng tin cậy, chẩn đoán bậc ba, kiểm định composition, tổng diện tích): 77/77 khớp sau làm tròn. Tất cả khớp. MPNc có AIC nhỏ nhất ở phân tích chính; MPNc ở phân tích chỉ có hệ số chặn.

## Dữ liệu

Nguồn: [AHMbook cswa.rda](https://github.com/cran/AHMbook/blob/8c006b16dae02d30a174c0dc754488c542e54b64/data/cswa.rda), commit `8c006b16dae02d30a174c0dc754488c542e54b64`; [tài liệu dữ liệu](https://github.com/cran/AHMbook/blob/8c006b16dae02d30a174c0dc754488c542e54b64/man/data_cswa.Rd). Bản gốc, tài liệu, SHA-256 và nguồn nằm trong `data/raw/`. Bộ giải mã dùng ở lần chạy này: `rda_reader`.

- 43 điểm; 128 lần khảo sát quan sát được; 42 điểm đủ 3 lần. Lần thiếu (`darty` lần 1) được loại khỏi likelihood, không mã hóa thành 0.
- Mảng nguồn có trục điểm × khoảng × lần; mã chuyển thành điểm × lần × khoảng.
- Diện tích được ghép theo tên điểm. Các điểm có thứ tự khác giữa bảng đếm và bảng covariate: `zeala`, `zealb`, `zeal4a`.
- Khoảng dài 2, 3, 5 phút; số phát hiện lần đầu theo khoảng [51, 14, 10], tổng 75 lượt phát hiện (không phải số cá thể riêng biệt).
- Tổng diện tích điểm 23.37 ha; exposure quan sát 69.68 ha-lần.
- Spot-mapping: 38 patch, 74 cá thể trên 66.69 ha = **1.109612 cá thể/ha** (chỉ dùng để mô tả, không vào likelihood).

## Phân tích chính với offset diện tích

Đặt $\mu_i=\rho A_i$. Với CRTE, cột mật độ là mật độ hiện diện trong một lần khảo sát $\nu_i/A_i$.

| Mô hình | logLik | AIC | ΔAIC | Trọng số | Mật độ/ha | r mỗi phút | c hoặc φ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MPNc | -173.833004 | 353.666008 | 0.000000 | 0.950507 | 1.083596 | 0.373167 | 0.456632 |
| CRTE | -176.788168 | 359.576336 | 5.910329 | 0.049493 | 1.088156 | 0.373191 | 0.424957 |
| MPN0 | -189.690004 | 383.380008 | 29.714000 | 0.000000 | 1.086521 | 0.373191 | — |
| MPN1 | -203.614175 | 411.228350 | 57.562343 | 0.000000 | 2.031182 | 0.086074 | — |

MPNc ước tính $\widehat\rho=1.083596$ cá thể/ha, chênh -2.34% so với spot-mapping; $\widehat c=0.456632$; xác suất phát hiện tích lũy 10 phút 0.990635. CRTE: mật độ hiện diện 1.088156, $\widehat\phi=0.424957$, mật độ superpopulation 2.560623 (CI log 1.684–3.895); đây là estimand khác. MPN1: mật độ 2.031182 (+83.1% so với spot-mapping), r = 0.086074.

## Khoảng tin cậy

CI chính được dựng trên thang log (mật độ) hoặc logit (xác suất) rồi biến đổi ngược, đúng như mô tả trong bản thảo. `confidence_intervals.csv` giữ thêm Wald trên thang gốc chỉ để đối chiếu.

| Mô hình | Tham số | Ước lượng | Hiệp phương sai | CI 95% (log/logit) |
| --- | --- | --- | --- | --- |
| MPN0 | rho | 1.087 | cluster | (0.773, 1.528) |
| MPN0 | r | 0.373 | cluster | (0.295, 0.458) |
| MPN1 | rho | 2.031 | model | (1.473, 2.801) |
| MPN1 | r | 0.086 | model | (0.062, 0.119) |
| MPNc | rho | 1.084 | model | (0.793, 1.481) |
| MPNc | r | 0.373 | model | (0.296, 0.457) |
| MPNc | c | 0.457 | model | (0.289, 0.635) |
| CRTE | rho | 1.088 | model | (0.801, 1.479) |
| CRTE | r | 0.373 | model | (0.296, 0.457) |
| CRTE | phi | 0.425 | model | (0.277, 0.587) |

MPN0 dùng sandwich gom cụm theo điểm; các mô hình khác dùng nghịch đảo observed information. Tỷ số SE sandwich/SE information của MPN0 là 1.504 cho mật độ và 1.016 cho r: phụ thuộc giữa các lần làm tăng chủ yếu phương sai của abundance, phù hợp với Corollary về hiệu chỉnh hạng một (đây là quan sát mô tả; thiết kế có offset nên công thức cân bằng không áp dụng chính xác).

## Kiểm định thành phần removal

Gộp qua các lần quan sát, số đếm theo khoảng là $\mathrm{Multinomial}(\sum d;\mathbf q(r))$ dưới mọi mô hình ứng viên. $\widehat r$ từ composition = 0.373191; kỳ vọng (45.96, 22.42, 6.62) so với quan sát [51, 14, 10]. Pearson $X^2=5.444$ (p = 0.0196), $G^2=5.681$ (p = 0.0171), 1 bậc tự do; p bootstrap tham số cho $G^2$ (4999 lần, seed 20261012) = 0.0222. Mô hình hazard hằng khớp kém với composition gộp, gợi ý dị biệt khả năng phát hiện.

## Chẩn đoán bậc ba

Tỷ số residual theo diện tích trên 42 điểm đủ ba lần: 0.332202 khi dùng ước lượng MPNc (dự đoán $\widehat c$ = 0.456632); 0.326303 khi dùng ước lượng CRTE (dự đoán $\widehat\phi^2$ = 0.180589). Bản gộp trừ trung bình theo lần (chỉ để audit): 0.376415. Đây là thống kê mô tả, không phải kiểm định.

## Phân tích độ nhạy: chỉ có hệ số chặn

Cột abundance có đơn vị cá thể/điểm/lần khảo sát.

| Mô hình | AIC | Trọng số | Số cá thể/điểm/lần | r mỗi phút | c hoặc φ |
| --- | --- | --- | --- | --- | --- |
| MPNc | 356.120210 | 0.836938 | 0.589969 | 0.371773 | 0.456554 |
| CRTE | 359.391472 | 0.163060 | 0.591729 | 0.373191 | 0.427994 |
| MPN0 | 381.743239 | 0.000002 | 0.591475 | 0.373191 | — |
| MPN1 | 410.844304 | 0.000000 | 1.100106 | 0.086601 | — |

## Kiểm tra số

50/50 lần tối ưu hội tụ; 8/8 phép khớp có CI. Chênh log-likelihood lớn nhất giữa các điểm khởi đầu 1.42e-13; |score| lớn nhất tại nghiệm 1.72e-06; thay đổi logLik khi siết ngưỡng cắt đuôi từ $10^{-12}$ xuống $10^{-14}$ và thêm 16 số hạng: 0. Đây là kiểm tra số, không phải chứng minh tối ưu toàn cục.

`results/paper_comparison.csv` liệt kê từng phép đối chiếu với dung sai làm tròn; `results/run_manifest.json` lưu phiên bản phần mềm, cấu hình và hash; `results/site_loglik.csv` lưu đóng góp theo điểm.

## Phạm vi

Các điểm được coi là độc lập dù thuộc 38 patch. Trọng số AIC không phải xác suất cơ chế sinh học đúng; dữ liệu không có định danh cá thể nên $c$ không xác nhận trực tiếp cá thể quay lại.

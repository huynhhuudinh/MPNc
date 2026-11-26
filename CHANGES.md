# Thay đổi so với gói gốc (2026-10-12)

## Kết luận kiểm tra

Không tìm thấy lỗi toán học trong các likelihood, thuật toán khớp, khoảng tin cậy,
chẩn đoán bậc ba hay phần tổng hợp mô phỏng. Các kiểm tra đã làm:

- Bốn likelihood được so với phép cộng trực tiếp trên mọi cách phân bổ multinomial
  gốc (kiểm thử sẵn có), và với phân phối bivariate Poisson khi có hai lần.
- Nghiệm profile `working_fit` trùng cực đại số của likelihood độc lập (sai khác
  log-likelihood < 1e-9).
- Bộ sinh dữ liệu cho đúng trung bình, hiệp phương sai và cross-cumulant bậc ba lý
  thuyết của cả MPNc và CRTE (400.000 site).
- Dữ liệu chuẩn hóa được giải mã lại độc lập từ `cswa.rda` (không qua `rdata`) và
  trùng khớp tuyệt đối, kể cả phép ghép ba điểm cuối theo tên.
- Application chạy lại: trùng lần chạy gốc đến ≥ 6 chữ số; 77/77 con số của Mục 7
  trong bản thảo khớp sau làm tròn.
- Mô phỏng chạy lại đúng ba ô đã công bố với cùng seed 20260910 (500 replicate mỗi
  ô): Bảng 3, tỷ lệ chọn đúng 58,3% (MCSE 0,0220), ΔAIC trung vị 0,172 và R̂₃
  (0,183; SD 0,300) khớp tuyệt đối; kết quả hai ô của Bảng 5 nằm trong
  `verification/`.

## Thay đổi

1. **`rda_reader.py` (mới).** Bộ giải mã RDX2/RDX3 không phụ thuộc gói ngoài.
   `run_application.read_source` dùng `rdata` nếu có, nếu không thì dùng bộ này;
   `run_manifest.json` ghi lại bộ giải mã đã dùng. Lý do: `rdata`/`xarray` không
   phải lúc nào cũng cài được, và cần một kiểm tra độc lập cho bước đọc dữ liệu.
2. **Kiểm định composition (`composition_gof`, mới).** Bản thảo đã sửa báo cáo
   X², G² và p cho composition gộp (Mục 7.3) nhưng code gốc không tính. Hàm mới ghi
   `results/composition_gof.json`, kèm p bootstrap tham số (4.999 lần, seed
   20261012) vì ô kỳ vọng nhỏ nhất chỉ 6,6.
3. **`build_report.py` viết lại.**
   - Bỏ mục "Điểm cần sửa: CI" đã lỗi thời (bản thảo hiện dùng CI log/logit).
   - Mọi câu trong báo cáo được tính từ kết quả thay vì viết cứng (diện tích,
     điểm ghép lại, chiều so sánh SE sandwich/information, mô hình có AIC nhỏ nhất).
   - Đối chiếu mở rộng từ 22 lên 77 con số, mỗi số có dung sai theo số chữ số
     được in; báo lỗi rõ nếu có số không khớp.
4. **`paper_simulation_outputs.py` (mới).** Phần Data Availability nói bảng và
   hình mô phỏng được tái tạo từ bản tóm tắt, nhưng gói gốc không có script cho
   hai hình `full_misspecification`, `full_model_selection` và các bảng gộp.
   Script mới dựng chúng từ thư mục kết quả của `softclosure.py`.
5. **Kiểm thử: 7 → 10.** Thêm: bộ giải mã độc lập so với dữ liệu chuẩn hóa;
   kiểm định composition (giá trị trên dữ liệu thật và trường hợp khớp gần hoàn
   hảo); profile working trùng MLE số.
6. **README.** Sửa tiêu đề theo bản thảo, bỏ đường dẫn máy cá nhân, mô tả tệp mới,
   cách dựng bảng/hình của bài và thời gian ước tính cho lần chạy đầy đủ.

Các tệp `softclosure.py` và `rda_reader.py`/`run_application.py` không thay đổi
bất kỳ công thức likelihood hay thiết lập tối ưu nào; kết quả số giữ nguyên.

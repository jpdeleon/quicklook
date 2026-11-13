# Changelog

## [1.10.1](https://github.com/jpdeleon/quicklook/compare/v1.10.0...v1.10.1) (2026-11-13)


### Performance Improvements

* cache TLS summary rows so /tls-summary stops hanging the GUI ([2d42a47](https://github.com/jpdeleon/quicklook/commit/2d42a47e32b0c573abcf6e69e30321965cfa1cf7))
* cache TLS summary rows so /tls-summary stops hanging the GUI ([ac45977](https://github.com/jpdeleon/quicklook/commit/ac4597755bdc9307acd3d3b011ef684adf2d8c16))

## [1.10.0](https://github.com/jpdeleon/quicklook/compare/v1.9.0...v1.10.0) (2026-11-12)


### Features

* add TOI-519 sectors 7 and 8 to the TGLC benchmark targets ([fcdff17](https://github.com/jpdeleon/quicklook/commit/fcdff175570c6b9b334165e9e167206e9d865a99))
* benchmark local TGLC light curves against the MAST HLSP ([57ba4f3](https://github.com/jpdeleon/quicklook/commit/57ba4f36c3d4f380aed4ed6c29e0282c22c30010))
* export TGLC quality flags and mask them in ql's local TGLC path ([6c783d7](https://github.com/jpdeleon/quicklook/commit/6c783d7a446938053b168ec415098e52d1d5cb00))
* expose the field-star prior in get_tglc_lc and the TGLC benchmark ([728f6f8](https://github.com/jpdeleon/quicklook/commit/728f6f803c17220740cb3e12c957e0684df2b4c3))


### Bug Fixes

* extract the requested target in get_tglc_lc, not the brightest star ([d1f24ac](https://github.com/jpdeleon/quicklook/commit/d1f24acd02c6ad7a7e171413ad3d63a7af93e7a1))

## [1.9.0](https://github.com/jpdeleon/quicklook/compare/v1.8.0...v1.9.0) (2026-11-11)


### Features

* extend TLS period ceiling to ~baseline and add catalogue epoch check ([242dd97](https://github.com/jpdeleon/quicklook/commit/242dd97e1424246aee03a804154163b53fb0ef36))
* extend TLS period ceiling to ~baseline and add catalogue epoch check ([c7f0db9](https://github.com/jpdeleon/quicklook/commit/c7f0db99de7a15e8f95b1953cec5bd37a7f10613))

## [1.8.0](https://github.com/jpdeleon/quicklook/compare/v1.7.0...v1.8.0) (2026-11-09)


### Features

* neighbor_check option; save Prot_ls and clear failed GLS ([1dc5b1f](https://github.com/jpdeleon/quicklook/commit/1dc5b1f76ccf296edd5f5053aa8707cb2a5c5f46))
* six detection/vetting improvements from gallery false-positive review ([7773344](https://github.com/jpdeleon/quicklook/commit/77733443dcdcd0e970ca416b7794f17bc5624355))


### Bug Fixes

* guard against missing per_transit_count in read-tls summary ([102b089](https://github.com/jpdeleon/quicklook/commit/102b08936b044a584b951cb894e2068a1ec9900a))
* penalize depth_ratio symmetrically in select_best_window ([933df68](https://github.com/jpdeleon/quicklook/commit/933df682350e426b0a8333b7b1cf69d25883e76a))

## [1.7.0](https://github.com/jpdeleon/quicklook/compare/v1.6.0...v1.7.0) (2026-11-05)


### Features

* add --sectors to stitch an explicit set of sectors ([ee3af86](https://github.com/jpdeleon/quicklook/commit/ee3af86038a0ba892a44245da1db14d317a2559b))
* add --sectors to stitch an explicit set of sectors ([508801a](https://github.com/jpdeleon/quicklook/commit/508801ab8025421daeb85e08ebf5719cb36c77ab))
* support --sector all in the run CLI ([dd853f1](https://github.com/jpdeleon/quicklook/commit/dd853f1fc81ad89636b30c5aa44615a0074d503b))
* support --sector all in the run CLI ([cb6ff3d](https://github.com/jpdeleon/quicklook/commit/cb6ff3df47da34e1e393552f203413faea52f2be))


### Bug Fixes

* do not abort run when Simbad query fails ([472b5f3](https://github.com/jpdeleon/quicklook/commit/472b5f378893ce603dc5055c9c936f7c85a12f07))
* do not abort run when Simbad query fails ([3948602](https://github.com/jpdeleon/quicklook/commit/39486023dfb5d201e35069603bcda45cf7655c49))

FROM genomic-aging-lfs:3.0-callers-permfix2 AS tests
USER root
COPY aging /opt/aging/aging
COPY resources /opt/aging/resources
COPY tests /opt/aging/tests
COPY prepare_cohort_config.py /opt/aging/prepare_cohort_config.py
COPY Dockerfile /opt/aging/Dockerfile
COPY finalize_bqr_feasibility.py /opt/aging/finalize_bqr_feasibility.py
COPY prepare_final_run.sh /opt/aging/prepare_final_run.sh
RUN chmod -R a+rX /opt/aging /opt/sage
USER mambauser
WORKDIR /opt/aging
RUN /usr/local/bin/_entrypoint.sh python -m unittest discover -s /opt/aging/tests -q

# The runnable image deliberately contains only the staged caller implementation.
# Legacy union/recount modules remain in the private source package for regression
# evidence but are not present in the production filesystem.
FROM genomic-aging-lfs:3.0-callers-permfix2
USER root
RUN mkdir -p /opt/aging/aging /opt/aging/resources
COPY --from=tests /opt/aging/aging/__init__.py /opt/aging/aging/common.py \
  /opt/aging/aging/contracts.py /opt/aging/aging/compact_qc.py \
  /opt/aging/aging/runtime.py /opt/aging/aging/stages.py \
  /opt/aging/aging/stage_contract.py /opt/aging/aging/stage_export.py \
  /opt/aging/aging/staged_production.py /opt/aging/aging/workbook_input.py \
  /opt/aging/aging/bqr_decision.py /opt/aging/aging/protocol.py \
  /opt/aging/aging/staged_inputs.py /opt/aging/aging/
COPY --from=tests /opt/aging/resources/compartments.GRCh38.json /opt/aging/resources/
RUN chmod -R a+rX /opt/aging /opt/sage
RUN /usr/local/bin/_entrypoint.sh python -c "import aging.staged_production, aging.stage_export" && \
  test ! -e /opt/aging/aging/reuse_callers.py && test ! -e /opt/aging/aging/single_pass.py
USER mambauser
WORKDIR /opt/aging
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "-m", "aging.staged_production"]

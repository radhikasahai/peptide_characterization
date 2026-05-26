// Production API URL is injected by .github/workflows/pages.yml from the
// PEPTIDE_API_URL repository variable. Local dev uses localhost:8000.
(function () {
  if (window.PEPTIDE_API_URL) return;

  var host = window.location.hostname;
  if (host === "localhost" || host === "127.0.0.1") {
    window.PEPTIDE_API_URL = "http://localhost:8000";
    return;
  }

  window.PEPTIDE_API_URL = "__PEPTIDE_API_URL__";
})();

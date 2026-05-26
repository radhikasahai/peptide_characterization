// Production API URL: GitHub Actions injects __PEPTIDE_API_URL__ from repo variable,
// or set via ?api=https://your-api.onrender.com or the in-page settings form.
(function () {
  if (window.PEPTIDE_API_URL) return;

  var params = new URLSearchParams(window.location.search);
  var fromQuery = params.get("api");
  if (fromQuery) {
    window.PEPTIDE_API_URL = fromQuery.replace(/\/$/, "");
    try {
      localStorage.setItem("PEPTIDE_API_URL", window.PEPTIDE_API_URL);
    } catch (e) {}
    return;
  }

  try {
    var stored = localStorage.getItem("PEPTIDE_API_URL");
    if (stored) {
      window.PEPTIDE_API_URL = stored.replace(/\/$/, "");
      return;
    }
  } catch (e) {}

  var host = window.location.hostname;
  if (host === "localhost" || host === "127.0.0.1") {
    window.PEPTIDE_API_URL = "http://localhost:8000";
    return;
  }

  window.PEPTIDE_API_URL = "__PEPTIDE_API_URL__";
})();

window.peptideApiConfigured = function () {
  return (
    window.PEPTIDE_API_URL &&
    window.PEPTIDE_API_URL.indexOf("__PEPTIDE_API_URL__") === -1
  );
};

(function () {
  function getCard(agent) {
    return document.querySelector(
      '[data-agent="' + agent.toLowerCase() + '"]'
    );
  }

  function updateAgentCard(agent, status, message) {
    const card = getCard(agent);

    if (!card) {
      return;
    }

    const normalizedStatus = String(status || "waiting").toLowerCase();
    const badge = card.querySelector(".agent-status");
    const messageElement = card.querySelector(".agent-message");

    card.classList.remove(
      "agent-waiting",
      "agent-running",
      "agent-completed",
      "agent-warning",
      "agent-failed"
    );

    card.classList.add("agent-" + normalizedStatus);

    if (badge) {
      badge.textContent =
        normalizedStatus.charAt(0).toUpperCase() +
        normalizedStatus.slice(1);
    }

    if (messageElement && message) {
      messageElement.textContent = message;
    }
  }

  function updatePipelineStatus(payload) {
    if (!payload || !payload.agent) {
      return;
    }

    updateAgentCard(
      payload.agent,
      payload.status,
      payload.message
    );
  }

  document.addEventListener("shiny:connected", function () {
    ["go", "kegg", "gsva", "chea"].forEach(function (agent) {
      updateAgentCard(
        agent,
        "waiting",
        "Waiting for analysis to begin."
      );
    });
  });

  if (window.Shiny) {
    Shiny.addCustomMessageHandler(
      "agent-status",
      updatePipelineStatus
    );
  }

  window.updateAgentCard = updateAgentCard;
})();

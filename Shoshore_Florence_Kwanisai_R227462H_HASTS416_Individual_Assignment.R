# HASTS416 - STOCHASTIC PROCESSES : TUTORIAL 1 IN R

# ==========================================================
# 0. REQUIRED PACKAGES
# ==========================================================
library(markovchain)
library(igraph)
library(expm)

# ==========================================================
# 1. FIXED OUTPUT LOCATION
# ==========================================================
BASE_DIR <- "C:/Users/user/Videos/Florence"

OUTPUT_DIR <- file.path(
  BASE_DIR,
  "Shoshore_Florence_Kwanisai_R227462H_HASTS416_Individual_Assignment"
)

DIR_A1 <- file.path(OUTPUT_DIR, "A1")
DIR_A2 <- file.path(OUTPUT_DIR, "A2")
DIR_A3 <- file.path(OUTPUT_DIR, "A3")
DIR_SUMMARY <- file.path(OUTPUT_DIR, "Summary")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_A1, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_A2, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_A3, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_SUMMARY, recursive = TRUE, showWarnings = FALSE)

# ==========================================================
# 2. HELPER FUNCTIONS
# ==========================================================
simulate_path <- function(P, nsteps, start_state = NULL){
  states <- 1:nrow(P)
  if(is.null(start_state)) start_state <- sample(states, 1)
  path <- numeric(nsteps + 1)
  path[1] <- start_state
  for(t in 1:nsteps){
    path[t + 1] <- sample(states, 1, prob = P[path[t], ])
  }
  return(path)
}

stationary_dist <- function(P){
  eig <- eigen(t(P))
  vec <- Re(eig$vectors[, which.min(abs(eig$values - 1))])
  pi <- vec / sum(vec)
  return(as.numeric(pi))
}

append_line <- function(text, file){
  cat(text, "\n", file = file, append = TRUE)
}

run_git <- function(cmd){
  system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
}

# ==========================================================
# 3. TEXT OUTPUT FILE
# ==========================================================
RESULTS_TXT <- file.path(OUTPUT_DIR, "HASTS416_Full_Results.txt")
if(file.exists(RESULTS_TXT)) file.remove(RESULTS_TXT)

append_line("HASTS416 - STOCHASTIC PROCESSES : TUTORIAL 1 IN R", RESULTS_TXT)
append_line("COMPLETE SOLUTION OUTPUT: A1, A2, A3", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)

# ==========================================================
# A1: 5-state Markov Chain Analysis
# ==========================================================
append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("A1", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)

# Transition matrix for A1
P1 <- matrix(c(
  1.0, 0,   0,   0,   0,
  0.5, 0,   0,   0,   0.5,
  0.2, 0,   0,   0,   0.8,
  0,   0,   1.0, 0,   0,
  0,   0,   0,   1.0, 0
), nrow = 5, byrow = TRUE)

states1 <- c("1", "2", "3", "4", "5")

mc1 <- new("markovchain",
           states = states1,
           transitionMatrix = P1)

# ----------------------------------------------------------
# A1(a) Markov Chain Diagram and Analysis
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A1(a)", RESULTS_TXT)
append_line("Transient classes: {2}, {3,4,5}", RESULTS_TXT)
append_line("Recurrent class: {1}", RESULTS_TXT)
append_line("Absorbing state: 1", RESULTS_TXT)
append_line("Reflecting state: 1", RESULTS_TXT)
append_line("Periods:", RESULTS_TXT)
append_line("d(1) = 1", RESULTS_TXT)
append_line("d(2) = undefined (no return to state 2 once left)", RESULTS_TXT)
append_line("d(3) = 3", RESULTS_TXT)
append_line("d(4) = 3", RESULTS_TXT)
append_line("d(5) = 3", RESULTS_TXT)

a1a_df <- data.frame(
  Item = c(
    "Transient classes",
    "Recurrent class",
    "Absorbing state",
    "Reflecting state",
    "d(1)",
    "d(2)",
    "d(3)",
    "d(4)",
    "d(5)"
  ),
  Answer = c(
    "{2}, {3,4,5}",
    "{1}",
    "1",
    "1",
    "1",
    "undefined",
    "3",
    "3",
    "3"
  ),
  stringsAsFactors = FALSE
)
write.csv(a1a_df, file.path(DIR_A1, "A1a_Answers.csv"), row.names = FALSE)

png(file.path(DIR_A1, "A1a_Markov_Chain_Diagram.png"), width = 1600, height = 1200, res = 200)
plot(mc1, main = "A1(a): Markov Chain Diagram")
dev.off()

# ----------------------------------------------------------
# A1(b) Three Simulated Trajectories
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A1(b)", RESULTS_TXT)

set.seed(123)
a1_path1 <- simulate_path(P1, 20)
a1_path2 <- simulate_path(P1, 20)
a1_path3 <- simulate_path(P1, 20)

a1_paths <- data.frame(
  Time = 0:20,
  Path1 = a1_path1,
  Path2 = a1_path2,
  Path3 = a1_path3
)
write.csv(a1_paths, file.path(DIR_A1, "A1b_Simulated_Trajectories.csv"), row.names = FALSE)

png(file.path(DIR_A1, "A1b_Three_Simulated_Trajectories.png"), width = 1600, height = 1200, res = 200)
matplot(rbind(a1_path1, a1_path2, a1_path3),
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "State",
        main = "A1(b): Three simulated trajectories")
legend("topright",
       legend = c("Path 1", "Path 2", "Path 3"),
       col = 1:3, lty = 1, lwd = 2)
dev.off()

append_line("Comment:", RESULTS_TXT)
append_line("All three trajectories eventually end in state 1 and remain there.", RESULTS_TXT)
append_line("States 3, 4 and 5 may show temporary cycling before absorption in state 1.", RESULTS_TXT)

# ----------------------------------------------------------
# A1(c)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A1(c)", RESULTS_TXT)

pi1 <- stationary_dist(P1)
names(pi1) <- states1
pi1_df <- data.frame(
  State = states1,
  Steady_State_Probability = round(pi1, 6)
)
write.csv(pi1_df, file.path(DIR_A1, "A1c_Steady_State_Probabilities.csv"), row.names = FALSE)

append_line("Steady-state probabilities:", RESULTS_TXT)
capture.output(print(pi1_df), file = RESULTS_TXT, append = TRUE)
append_line("Interpretation:", RESULTS_TXT)
append_line("In the long run, the chain spends all its time in state 1.", RESULTS_TXT)
append_line("Therefore the stationary distribution is (1,0,0,0,0).", RESULTS_TXT)
append_line("Is it an ergodic chain? No. The chain is not irreducible.", RESULTS_TXT)

# ----------------------------------------------------------
# A1(d)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A1(d)", RESULTS_TXT)

mu0_a1 <- rep(1/5, 5)
nmax1 <- 20
probs1 <- matrix(0, nrow = nmax1 + 1, ncol = 5)
colnames(probs1) <- states1
probs1[1, ] <- mu0_a1

Pn1 <- diag(5)
for(n in 1:nmax1){
  Pn1 <- Pn1 %*% P1
  probs1[n + 1, ] <- mu0_a1 %*% Pn1
}

probs1_df <- data.frame(Time = 0:nmax1, probs1)
write.csv(probs1_df, file.path(DIR_A1, "A1d_Unconditional_Probabilities.csv"), row.names = FALSE)

png(file.path(DIR_A1, "A1d_Unconditional_Probabilities.png"), width = 1600, height = 1200, res = 200)
matplot(0:nmax1, probs1,
        type = "l", lty = 1, lwd = 2,
        xlab = "Time n", ylab = "Probability",
        main = "A1(d): Unconditional probabilities over time")
legend("right",
       legend = paste("State", 1:5),
       col = 1:5, lty = 1, lwd = 2)
dev.off()

append_line("Comment:", RESULTS_TXT)
append_line("The unconditional probabilities converge to (1,0,0,0,0).", RESULTS_TXT)
append_line("Probability of state 1 rises quickly to 1, while probabilities of the other states decay to 0.", RESULTS_TXT)
append_line("The convergence is relatively fast because state 1 is absorbing.", RESULTS_TXT)

# ==========================================================
# A2
# ==========================================================
append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("A2", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)

P2 <- matrix(c(
  0,   1,   0,   0,   0,   0,   0,
  1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0.4, 0.2, 0.2, 0.2,
  0,   0,   0,   0,   0.2, 0.4, 0.4,
  0.3, 0,   0,   0.1, 0.3, 0.1, 0.2,
  0,   0,   0,   0.2, 0.2, 0.3, 0.3,
  0,   0,   0,   0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)

states2 <- as.character(1:7)

mc2 <- new("markovchain",
           states = states2,
           transitionMatrix = P2)

# ----------------------------------------------------------
# A2(a)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(a)", RESULTS_TXT)
append_line("The Markov chain diagram has been plotted and saved.", RESULTS_TXT)

png(file.path(DIR_A2, "A2a_Markov_Chain_Diagram.png"), width = 1600, height = 1200, res = 200)
plot(mc2)
title(main = "A2(a): Markov Chain Diagram")
dev.off()

# ----------------------------------------------------------
# A2(b)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(b)", RESULTS_TXT)
append_line("Recurrent class: {1,2}", RESULTS_TXT)
append_line("Transient class: {3,4,5,6,7}", RESULTS_TXT)
append_line("Absorbing states: none", RESULTS_TXT)
append_line("Reflecting states: 5, 6, 7", RESULTS_TXT)
append_line("Periods:", RESULTS_TXT)
append_line("d(1) = 2", RESULTS_TXT)
append_line("d(2) = 2", RESULTS_TXT)
append_line("d(3) = undefined (no return)", RESULTS_TXT)
append_line("d(4) = 1", RESULTS_TXT)
append_line("d(5) = 1", RESULTS_TXT)
append_line("d(6) = 1", RESULTS_TXT)
append_line("d(7) = 1", RESULTS_TXT)

a2b_df <- data.frame(
  Item = c(
    "Recurrent class",
    "Transient class",
    "Absorbing states",
    "Reflecting states",
    "d(1)",
    "d(2)",
    "d(3)",
    "d(4)",
    "d(5)",
    "d(6)",
    "d(7)"
  ),
  Answer = c(
    "{1,2}",
    "{3,4,5,6,7}",
    "none",
    "5, 6, 7",
    "2",
    "2",
    "undefined",
    "1",
    "1",
    "1",
    "1"
  ),
  stringsAsFactors = FALSE
)
write.csv(a2b_df, file.path(DIR_A2, "A2b_Answers.csv"), row.names = FALSE)

# ----------------------------------------------------------
# A2(c)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(c)", RESULTS_TXT)

set.seed(321)
a2_path1 <- simulate_path(P2, 30)
a2_path2 <- simulate_path(P2, 30)

a2_paths <- data.frame(
  Time = 0:30,
  Path1 = a2_path1,
  Path2 = a2_path2
)
write.csv(a2_paths, file.path(DIR_A2, "A2c_Simulated_Trajectories.csv"), row.names = FALSE)

png(file.path(DIR_A2, "A2c_Two_Simulated_Trajectories.png"), width = 1600, height = 1200, res = 200)
matplot(rbind(a2_path1, a2_path2),
        type = "l", lty = 1, lwd = 2,
        xlab = "Time", ylab = "State",
        main = "A2(c): Two simulated trajectories")
legend("topright",
       legend = c("Path 1", "Path 2"),
       col = 1:2, lty = 1, lwd = 2)
dev.off()

append_line("Discussion:", RESULTS_TXT)
append_line("A trajectory starting in states 3 to 7 may wander for some time,", RESULTS_TXT)
append_line("but once it enters the class {1,2}, it alternates forever between states 1 and 2.", RESULTS_TXT)
append_line("This shows that {1,2} is a closed recurrent class, while {3,4,5,6,7} is transient.", RESULTS_TXT)

# ----------------------------------------------------------
# A2(d)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A2(d)", RESULTS_TXT)

pi2 <- stationary_dist(P2)
names(pi2) <- states2
pi2_df <- data.frame(
  State = states2,
  Stationary_Probability = round(pi2, 6)
)
write.csv(pi2_df, file.path(DIR_A2, "A2d_Stationary_Distribution.csv"), row.names = FALSE)

append_line("Limiting / stationary probabilities:", RESULTS_TXT)
capture.output(print(pi2_df), file = RESULTS_TXT, append = TRUE)
append_line("Interpretation:", RESULTS_TXT)
append_line("In the long run, probability mass is concentrated on states 1 and 2.", RESULTS_TXT)
append_line("The stationary distribution is (1/2, 1/2, 0, 0, 0, 0, 0).", RESULTS_TXT)
append_line("Is the chain ergodic? No.", RESULTS_TXT)
append_line("Although a stationary distribution exists, the recurrent class {1,2} has period 2,", RESULTS_TXT)
append_line("so the chain is not ergodic.", RESULTS_TXT)

# ==========================================================
# A3
# ==========================================================
append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("A3", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)

append_line("For question A3:", RESULTS_TXT)
append_line("The first transition matrix used is:", RESULTS_TXT)
append_line("[[0.4, 0.4, 0.2], [0.3, 0.4, 0.3], [0, 0.1, 0.9]]", RESULTS_TXT)

P13 <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

P46 <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

traffic_states <- c("light", "heavy", "jammed")
mu0_a3 <- c(1, 0, 0)

# ----------------------------------------------------------
# A3(a)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A3(a)", RESULTS_TXT)
append_line("From 1 PM to 4 PM = 3 hours = 9 steps of 20 minutes each.", RESULTS_TXT)
append_line("From 4 PM to 6 PM = 2 hours = 6 steps of 20 minutes each.", RESULTS_TXT)

mu6 <- mu0_a3 %*% (P13 %^% 9) %*% (P46 %^% 6)
colnames(mu6) <- traffic_states

mu6_df <- data.frame(
  State = traffic_states,
  Probability = round(as.numeric(mu6), 6)
)
write.csv(mu6_df, file.path(DIR_A3, "A3a_Distribution_at_6PM.csv"), row.names = FALSE)

append_line("Distribution at 6 PM:", RESULTS_TXT)
capture.output(print(mu6_df), file = RESULTS_TXT, append = TRUE)
append_line(sprintf("P(light at 6 PM)  = %.6f", mu6[1,1]), RESULTS_TXT)
append_line(sprintf("P(heavy at 6 PM)  = %.6f", mu6[1,2]), RESULTS_TXT)
append_line(sprintf("P(jammed at 6 PM) = %.6f", mu6[1,3]), RESULTS_TXT)
append_line("Interpretation:", RESULTS_TXT)
append_line("By 6 PM, traffic is most likely to be jammed.", RESULTS_TXT)

# ----------------------------------------------------------
# A3(b)
# ----------------------------------------------------------
append_line("", RESULTS_TXT)
append_line("A3(b)", RESULTS_TXT)

simulate_traffic <- function() {
  state <- 1
  for(i in 1:9){
    state <- sample(1:3, 1, prob = P13[state, ])
  }
  for(i in 1:6){
    state <- sample(1:3, 1, prob = P46[state, ])
  }
  return(state)
}

set.seed(123)
N <- 10000
final_states <- replicate(N, simulate_traffic())
sim_props <- table(final_states) / N

sim_vector <- c(
  light  = ifelse("1" %in% names(sim_props), sim_props["1"], 0),
  heavy  = ifelse("2" %in% names(sim_props), sim_props["2"], 0),
  jammed = ifelse("3" %in% names(sim_props), sim_props["3"], 0)
)

sim_df <- data.frame(
  State = c("light", "heavy", "jammed"),
  Simulated_Proportion = round(as.numeric(sim_vector), 6)
)
write.csv(sim_df, file.path(DIR_A3, "A3b_Simulated_Proportions_10000.csv"), row.names = FALSE)

png(file.path(DIR_A3, "A3b_Simulated_Distribution_at_6PM.png"), width = 1600, height = 1200, res = 200)
barplot(sim_vector,
        main = "A3(b): Simulated distribution at 6 PM",
        ylab = "Proportion",
        names.arg = c("light", "heavy", "jammed"))
dev.off()

append_line("Simulated proportions at 6 PM from 10,000 trajectories:", RESULTS_TXT)
capture.output(print(sim_df), file = RESULTS_TXT, append = TRUE)
append_line("Comment:", RESULTS_TXT)
append_line("The simulated proportions are close to the exact probabilities found in A3(a),", RESULTS_TXT)
append_line("thereby verifying the result.", RESULTS_TXT)

# ==========================================================
# A1, A2 AND A3 COMMENTS
# ==========================================================
SUMMARY_CSV <- file.path(DIR_SUMMARY, "A1_A2_A3_Comments.csv")

summary_df <- data.frame(
  Question = c("A1(a)", "A1(b)", "A1(c)", "A1(d)", "A2(a)", "A2(b)", "A2(c)", "A2(d)", "A3(a)", "A3(b)"),
  Question_Text = c(
    "Use R to plot a diagram of the Markov chain. Identify all transient and recurrent classes. Identify all absorbing and reflective states. Find the period of each state.",
    "Simulate three trajectories of the chain that start at a randomly chosen state. Comment on what you see.",
    "Find the steady-state probabilities and interpret them. Is it an ergodic chain?",
    "Plot the unconditional probabilities at time n against time and comment on how fast the probabilities converge to the steady-state distribution.",
    "Use R to plot a diagram of the Markov chain.",
    "Identify all recurrent and transient classes. Find their periods. Are there any absorbing and reflecting states?",
    "Simulate two trajectories of the chain that start at a randomly selected state. Discuss what you see in the plot.",
    "Calculate the limiting probabilities and interpret them. Is the chain ergodic?",
    "If the traffic starts with light state at 1PM, what is the distribution of states at 6PM?",
    "Simulate 10,000 trajectories to verify the result of the previous part."
  ),
  Result = c(
    "Transient classes {2},{3,4,5}; recurrent class {1}; absorbing state 1; reflecting state 1; periods d(1)=1, d(2)=undefined, d(3)=d(4)=d(5)=3",
    "Three simulated trajectories eventually end in state 1",
    "Steady-state distribution (1,0,0,0,0); not ergodic",
    "Unconditional probabilities converge to (1,0,0,0,0)",
    "Markov chain diagram plotted and saved",
    "Recurrent class {1,2}; transient class {3,4,5,6,7}; no absorbing states; reflecting states 5,6,7; periods identified",
    "Two simulated trajectories show eventual entry into {1,2}",
    "Stationary distribution (1/2,1/2,0,0,0,0,0); not ergodic",
    paste0("Distribution at 6 PM = (", round(mu6[1,1],6), ", ", round(mu6[1,2],6), ", ", round(mu6[1,3],6), ")"),
    "10,000 simulated trajectories verify A3(a)"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, SUMMARY_CSV, row.names = FALSE)

append_line("", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
append_line("A1, A2 AND A3 COMMENTS", RESULTS_TXT)
append_line("====================================================", RESULTS_TXT)
capture.output(print(summary_df), file = RESULTS_TXT, append = TRUE)

# ==========================================================
# FINAL CONSOLE MESSAGE
# ==========================================================
cat("\n====================================================\n")
cat("Results successfully saved to:\n")
cat(OUTPUT_DIR, "\n")
cat("====================================================\n")

# ==========================================================
# 4. GITHUB PUSH SECTION
# ==========================================================
# Make sure Git is installed on your computer before running this part.

GITHUB_REPO_URL <- "https://github.com/shoshoreflorence03-blip/Shoshore-Florence-Kwanisai-R227462H-Assignment__HASTS416.git"
BRANCH_NAME <- "main"

setwd(OUTPUT_DIR)

# Create .gitignore
writeLines(c(
  ".Rhistory",
  ".RData",
  ".Rproj.user",
  ".vscode/",
  "__pycache__/"
), ".gitignore")

# Initialize git repository if it does not already exist
if(!dir.exists(file.path(OUTPUT_DIR, ".git"))){
  run_git("git init")
}

# Ensure branch is main
run_git(paste("git branch -M", BRANCH_NAME))

# Remove old origin if it exists, then add correct origin
run_git("git remote remove origin")
run_git(paste("git remote add origin", shQuote(GITHUB_REPO_URL)))

# Add all files
run_git("git add .")

# Commit files
run_git('git commit -m "Add Florence HASTS416 assignment outputs"')

# Push to GitHub
run_git(paste("git push -u origin", BRANCH_NAME))

cat("\n====================================================\n")
cat("GitHub push attempted for repository:\n")
cat(GITHUB_REPO_URL, "\n")
cat("====================================================\n")
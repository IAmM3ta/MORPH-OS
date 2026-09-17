/** Morphogen-native ESC 1.0 L0 mapper for the 20 Hz Sync pump. */

export type GrayScottPreset =
  | 'mitosis' | 'solitons' | 'pulsing' | 'holes' | 'mazes'
  | 'fingerprint' | 'spirals' | 'worms' | 'coral' | 'uskate';

export type ResidualMode = 'L0' | 'L1' | 'L2' | 'L3' | 'auto';

export type EscPacket = {
  esc_version: '1.0';
  residual_mode?: ResidualMode;
  ontology_pointer: { class: 'GrayScott.TuringField'; preset: GrayScottPreset; prior_hash: string };
  manifold_anchor: { feed: number; kill: number; du: number; dv: number; seed: number; size: number };
  extropic_residual: {
    N_true: number; N_proxy_lei: number; C_score: number; L_score: number; Ex: number;
    mutations: { type: 'localized_defect'; cx: number; cy: number; magnitude: number; radius: number }[];
    homology?: { kind: 'superlevel_H0_periodic'; working_t: number; beta0: number; barcode: { birth: number; death: number; persist: number }[] };
  };
  temporal_delta: { time_step: number; delta_w1: number | null };
};

const PRESET_TABLE: Record<GrayScottPreset, { feed: number; kill: number; du: number; dv: number }> = {
  mitosis: { feed: 0.0367, kill: 0.0649, du: 0.16, dv: 0.08 },
  solitons: { feed: 0.0353, kill: 0.0653, du: 0.16, dv: 0.08 },
  pulsing: { feed: 0.025, kill: 0.06, du: 0.14, dv: 0.07 },
  holes: { feed: 0.039, kill: 0.058, du: 0.16, dv: 0.08 },
  mazes: { feed: 0.029, kill: 0.057, du: 0.16, dv: 0.08 },
  fingerprint: { feed: 0.026, kill: 0.061, du: 0.16, dv: 0.08 },
  spirals: { feed: 0.018, kill: 0.051, du: 0.16, dv: 0.08 },
  worms: { feed: 0.046, kill: 0.063, du: 0.16, dv: 0.08 },
  coral: { feed: 0.0545, kill: 0.062, du: 0.16, dv: 0.08 },
  uskate: { feed: 0.062, kill: 0.0609, du: 0.16, dv: 0.08 },
};

export function presetFromParams(feed: number, kill: number, eps = 1e-4): GrayScottPreset | null {
  let best: GrayScottPreset | null = null;
  let bestD = Infinity;
  (Object.keys(PRESET_TABLE) as GrayScottPreset[]).forEach((id) => {
    const p = PRESET_TABLE[id];
    const d = Math.abs(p.feed - feed) + Math.abs(p.kill - kill);
    if (d < bestD) { bestD = d; best = id; }
  });
  return bestD <= eps ? best : null;
}
